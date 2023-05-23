using Asap, AsapToolkit;
using Zygote
using kjlMakie; set_theme!(kjl_dark)

### Create a spaceframe

#meta parameters
begin
    nx = 15
    dx = 1500.
    ny = 35
    dy = 1500.
    dz = 3200.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
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

    hidedecorations!(ax); hidespines!(ax)

    e_init = linesegments!(e0,
        color = axf,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 3)

    fig
end

# parameters
positions = Asap.nodePositions(model)
Xo = positions[:, 1]; Yo = positions[:, 2]; Zo = positions[:, 3]
E = Steel_Nmm.E; A = sec.A

# collect parameters
p = TrussOptParams(model)

# objective function
function compliance(dX::Vector{Float64}, dY::Vector{Float64}, dZ::Vector{Float64}, p::TrussOptParams)

    X = Xo .+ dX
    Y = Yo .+ dY
    Z = Zo .+ dZ

    ks = [kglobal(X, Y, Z, E, A, id) for id in p.nodeids]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    U' * p.P[p.freeids]
end

# initial Δ
pzero = zero(Xo)

# test
@time compliance(pzero, pzero, pzero, p)

# closure function for Z only
Zcompliance(dZ::Vector{Float64}, p::TrussOptParams) = compliance(pzero, pzero, dZ, p)

# test
@time c = Zcompliance(pzero, p);

# gradient
@time dcdz = Zygote.gradient(var -> Zcompliance(var, p), pzero)[1];

# initial position and bounds
z0 = zero(Zo)
lb = -1500.
ub = 2500.

# assemble
func = Optimization.OptimizationFunction(Zcompliance, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(func, z0, p;
    lb = lb,
    ub = ub)

# storage/callback
begin
    losstrace = Vector{Float64}()
    ztrace = Vector{Vector{Float64}}()

    function cb(vars::Vector{Float64}, loss::Float64)
        push!(losstrace, loss)
        push!(ztrace, vars)
        false
    end
end

# solve
@time sol = Optimization.solve(prob,
    OptimizationNLopt.NLopt.LD_LBFGS(),
    reltol = 1e-3,
    callback = cb);

# assemble into new model
begin
    model2 = deepcopy(model)
    for (i, node) in enumerate(model2.nodes)
        node.position[3] += sol.u[i]
    end

    solve!(model2; reprocess = true)


    p1 = Point3.(getproperty.(model2.nodes, :position))
    e1 = vcat([p1[id] for id in getproperty.(model2.elements, :nodeIDs)]...)
    axf1 = getindex.(getproperty.(model2.elements, :forces), 2)
    cr1 = @lift(maximum(abs.(axf1)) .* (-1, 1) .* $crfac)
end

# compare
begin
    fig = Figure()
    ax = Axis3(fig[1,1],
        title = "Compliance: $(round(model.compliance))",
        aspect = :data)

    hidedecorations!(ax); hidespines!(ax)
    ax.titlevisible = true

    e_init = linesegments!(e0,
        color = axf,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 3)


    ax2 = Axis3(fig[1,2],
        title = "Compliance: $(round(model2.compliance))",
        aspect = :data)

    hidedecorations!(ax2); hidespines!(ax2)
    ax2.titlevisible = true

    e_sol = linesegments!(e1,
        color = axf1,
        colorrange = cr1,
        colormap = pink2blue,
        linewidth = 5)

    fig
end

### using partial variables
p = TrussOptParams(model)

iActive = [i for i in sf.ibottom if !in(i, sf.isupport)]
iActive = vec(sf.itop)

iActive = [i for i in collect(1:length(Zo)) if !in(i, sf.isupport)]

iActive = vec(sf.itop)
z0 = Zo[iActive]

function Zfreecompliance(Z::Vector{Float64}, p::TrussOptParams)
    Znew = updatevalues(Zo, iActive, Z)

    ks = [AsapToolkit.kglobal(Xo, Yo, Znew, E, A, id) for id in p.nodeids]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    U' * p.P[p.freeids]
end

@time Zfreecompliance(z0, p)
g = Zygote.gradient(var -> Zfreecompliance(var, p), z0)[1]

lb = z0 .- 1250
ub = z0 .+ 6000

func = Optimization.OptimizationFunction(Zfreecompliance, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(func, z0, p;
    lb = lb,
    ub = ub)

# storage/callback
begin
    losstrace = Vector{Float64}()
    valtrace = Vector{Vector{Float64}}()

    function cb(vars::Vector{Float64}, loss::Float64)
        push!(losstrace, loss)
        push!(valtrace, vars)
        false
    end
end

# solve
@time sol = Optimization.solve(prob,
    OptimizationNLopt.NLopt.LD_LBFGS(),
    reltol = 1e-3,
    callback = cb);

