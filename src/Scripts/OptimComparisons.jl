using Asap, AsapToolkit
using Zygote, LinearAlgebra, FiniteDiff
using kjlMakie; set_theme!(kjl_dark)

### Create a spaceframe

#meta parameters
begin
    nx = 10
    dx = 750.
    ny = 10
    dy = 1500.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :corner);
model = sf.truss;

#plot assets
begin
    crfac = Observable(0.5)
    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)
    axf = getindex.(getproperty.(model.elements, :forces), 2)
    cr = @lift(maximum(abs.(axf)) .* (-1, 1) .* $crfac)
end

#plot initial spaceframe
begin
    fig = Figure()
    ax = Axis3(fig[1,1],
        aspect = :data)

    gridtoggle!(ax)

    e_init = linesegments!(e0,
        color = axf,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 3)

    fig
end

#params
p = TrussOptParams(model)
positions = Asap.nodePositions(model)
Xo = positions[:, 1]; Yo = positions[:, 2]; Zo = positions[:, 3]
E = Steel_Nmm.E; A = sec.A

itopcorners = [sf.itop[1,1], sf.itop[1,end], sf.itop[end,1], sf.itop[end,end]]
iActive = [i for i in eachindex(model.nodes) if !in(i, sf.isupport)]

#initial variables
z0 = Zo[iActive]

#obj func
function Zfreecompliance(Z::Vector{Float64}, p::TrussOptParams)
    Znew = updatevalues(Zo, iActive, Z)

    ks = [AsapToolkit.kglobal(Xo, Yo, Znew, E, A, id) for id in p.nodeids]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    U' * p.P[p.freeids]
end

#obj func without adjoint
function ZfreecomplianceNoAdjoint(Z::Vector{Float64}, p::TrussOptParams)
    Znew = updatevalues(Zo, iActive, Z)

    ks = [AsapToolkit.kglobal(Xo, Yo, Znew, E, A, id) for id in p.nodeids]

    K = assembleglobalK(ks, p)

    U = K[p.freeids, p.freeids] \ p.P[p.freeids]

    U' * p.P[p.freeids]
end

@time z1 = Zfreecompliance(z0, p)
@time z1b = ZfreecomplianceNoAdjoint(z0, p)

@time dz1 = Zygote.gradient(var -> Zfreecompliance(var, p), z0)[1];
@time dz1b = Zygote.gradient(var -> ZfreecomplianceNoAdjoint(var, p), z0)[1];

#bounds
lb = z0 .- 2000
ub = z0 .+ 2000

newitops = [findfirst(i .== iActive) for i in itopcorners]
lb[newitops] .= z0[newitops]

#solve with adjoint
begin
    func1 = Optimization.OptimizationFunction(Zfreecompliance, Optimization.AutoZygote())
    prob1 = Optimization.OptimizationProblem(func1, z0, p;
        lb = lb,
        ub = ub)

    trace1 = Vector{Float64}()
    vals1 = Vector{Vector{Float64}}()

    function cb1(vars::Vector{Float64}, loss::Float64)
        push!(trace1, loss)
        push!(vals1, vars)
        false
    end
    @time sol1 = Optimization.solve(prob1,
        OptimizationNLopt.NLopt.LD_LBFGS(),
        reltol = 1e-3,
        callback = cb1);
end

#solve without adjoint
begin
    func2 = Optimization.OptimizationFunction(ZfreecomplianceNoAdjoint, Optimization.AutoZygote())
    prob2 = Optimization.OptimizationProblem(func2, z0, p;
        lb = lb,
        ub = ub)

    trace2 = Vector{Float64}()
    vals2 = Vector{Vector{Float64}}()

    function cb2(vars::Vector{Float64}, loss::Float64)
        push!(trace2, loss)
        push!(vals2, vars)
        false
    end
    @time sol2 = Optimization.solve(prob2,
        OptimizationNLopt.NLopt.LD_LBFGS(),
        reltol = 1e-3,
        callback = cb2);
end

#solve with finite difference
begin
    func3 = Optimization.OptimizationFunction(Zfreecompliance, Optimization.AutoFiniteDiff())
    prob3 = Optimization.OptimizationProblem(func3, z0, p;
        lb = lb,
        ub = ub)

    trace3 = Vector{Float64}()
    vals3 = Vector{Vector{Float64}}()

    function cb3(vars::Vector{Float64}, loss::Float64)
        push!(trace3, loss)
        push!(vals3, vars)
        false
    end
    @time sol3 = Optimization.solve(prob3,
        OptimizationNLopt.NLopt.LD_LBFGS(),
        reltol = 1e-3,
        callback = cb3);
end

#solve with gradient-free method difference
begin
    func4 = Optimization.OptimizationFunction(Zfreecompliance, Optimization.AutoFiniteDiff())
    prob4 = Optimization.OptimizationProblem(func4, z0, p;
        lb = lb,
        ub = ub)

    trace4 = Vector{Float64}()
    vals4 = Vector{Vector{Float64}}()

    function cb4(vars::Vector{Float64}, loss::Float64)
        push!(trace4, loss)
        push!(vals4, vars)
        false
    end
    @time sol4 = Optimization.solve(prob4,
        OptimizationNLopt.NLopt.LN_BOBYQA(),
        reltol = 1e-3,
        callback = cb4);
end

#convert to models
begin
    model1 = deepcopy(model)
    for (i, node) in enumerate(model1.nodes[iActive])
        node.position[3] = sol1.u[i]
    end

    solve!(model1; reprocess = true)

    model2 = deepcopy(model)
    for (i, node) in enumerate(model2.nodes[iActive])
        node.position[3] = sol2.u[i]
    end

    solve!(model2; reprocess = true)

    model3 = deepcopy(model)
    for (i, node) in enumerate(model3.nodes[iActive])
        node.position[3] = sol3.u[i]
    end

    solve!(model3; reprocess = true)

    model4 = deepcopy(model)
    for (i, node) in enumerate(model4.nodes[iActive])
        node.position[3] = sol4.u[i]
    end

    solve!(model4; reprocess = true)
end

#extract visualization assets
begin
    pts1 = Point3.(getproperty.(model1.nodes, :position))
    e1 = vcat([pts1[id] for id in p.nodeids]...)
    f1 = getindex.(getproperty.(model1.elements, :forces), 2)

    pts2 = Point3.(getproperty.(model2.nodes, :position))
    e2 = vcat([pts2[id] for id in p.nodeids]...)
    f2 = getindex.(getproperty.(model2.elements, :forces), 2)

    pts3 = Point3.(getproperty.(model3.nodes, :position))
    e3 = vcat([pts3[id] for id in p.nodeids]...)
    f3 = getindex.(getproperty.(model3.elements, :forces), 2)

    pts4 = Point3.(getproperty.(model4.nodes, :position))
    e4 = vcat([pts4[id] for id in p.nodeids]...)
    f4 = getindex.(getproperty.(model4.elements, :forces), 2)
end

begin
    fig = Figure()

    ax0 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0,
        color = axf,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 3)

    ax1 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax1); hidespines!(ax1)

    linesegments!(e1,
        # color = f1,
        # colorrange = cr,
        # colormap = pink2blue,
        color = blue,
        linewidth = 3)


    ax3 = Axis3(fig[1,3],
        aspect = :data)

    hidedecorations!(ax3); hidespines!(ax3)

    linesegments!(e3,
        # color = f3,
        # colorrange = cr,
        # colormap = pink2blue,
        color = :white,
        linewidth = 3)

    ax4 = Axis3(fig[1,4],
        aspect = :data)

    hidedecorations!(ax4); hidespines!(ax4)

    linesegments!(e4,
        # color = f4,
        # colorrange = cr,
        # colormap = pink2blue,
        color = green,
        linewidth = 3)


    axloss = Axis(fig[2,1:3],
        aspect = nothing,
        xlabel = "Iteration",
        ylabel = "Compliance [Nmm]")

    t1 = lines!(trace1,
        label = "AD + Adjoint.",
        linewidth = 3,
        color = blue)

    t2 = lines!(trace2,
        label = "AD Only",
        linewidth = 3,
        color = :lightblue)

    t3 = lines!(trace3,
        label = "Finite Diff.",
        linewidth = 3,
        color = :white)

    t4 = lines!(trace4,
        label = "Gradient Free",
        linewidth = 3,
        color = green)

    leg = axislegend(axloss,
        position = :rb)

    axtime = Axis(fig[2,4],
        aspect = nothing,
        ylabel = "Solve time [s]")

    times = barplot!(1:4, [sol1.solve_time, sol2.solve_time, sol3.solve_time, sol4.solve_time],
        color = [blue, :lightblue, :white, green])

    colsize!(fig.layout, 4, Aspect(2, 1.))

    fig
end