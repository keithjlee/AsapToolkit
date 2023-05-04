using DifferentialEquations

#differential eq. solution
p = rand(model.elements[:primary])
uends = uImag[[2, 6, 8, 12]]
xfine = range(0, p.length, 200)
ushapefunc = vcat([Asap.N(x, p.length) * uends for x in xfine]...)

klocal = Asap.localK(p)
uImag = [[0., 100., 0., 0., 0., deg2rad(20)]; zeros(6)]
fImag = klocal * uImag

v, m = fImag[[2, 6]]
EI = p.section.E * p.section.Izz

function ddu(u′, u, p, t)
    1 / EI * (v * t - m)
end

ndiscrete = 200
θₒ = uImag[6]
uₒ = uImag[2]
tspan = (0., p.length)
prob = SecondOrderODEProblem(ddu, θₒ, uₒ, tspan)

xset = Vector{Vector{Float64}}()
uset = Vector{Vector{Float64}}()
nset = Vector{Int64}()

for ndiscrete = 5:200

    push!(nset, ndiscrete)
    @time sol = DifferentialEquations.solve(prob, DPRKN6(); tstops = range(0, p.length, ndiscrete))

    u = getindex.(sol.u, 2)
    xrange = range(tspan..., length(u))

    push!(uset, u)
    push!(xset, collect(xrange))

end

i = Observable(1)


nsol = @lift("ODE solution at " * string(nset[$i]) * " discretizations")
xusol = @lift(Point2.(xset[$i], uset[$i]))

begin
    lw = 5
    fig = Figure(backgroundcolor = :black)
    ax = Axis(fig[1,1], aspect = 1,
        xlabel = "local x [mm]",
        title = nsol,
        ylabel = "Δ [mm]")

    lines!(xfine, ushapefunc,
        linewidth = lw,
        color = blue,
        label = "Hermitian shape function 3rd order")

    lines!(xusol,
            linewidth = lw,
            color = :white,
            label = "ODE Sol.")

    axislegend(ax, position = :lb)
    fig
end

for k = 1:length(nset)
    i[] = k
    sleep(.1)
end