using DifferentialEquations

#differential eq. solution
begin
    n1 = Node([0., 0., 0.], :pinned)
    n2 = Node([8000., 0., 0.], :fixed)

    p = Element(n1, n2, joist)
    # p.Ψ = 0.

    l = LineLoad(p, [0., 20., -20.])

    nodes = [n1, n2]; elements = [p]; loads = [l]

    model = Model(nodes, elements, loads)
    solve!(model)

    n = 20
    xlocal = first(p.LCS)
    init = p.nodeStart.position

    u = AsapToolkit.dofdisplacement(p; n = n)
    inc = range(0, p.length, n)
    df = Observable(10.)

    xyz = @lift(Point3.(eachcol(hcat([init .+ xlocal .* i .- $df .* sum(disp .* p.LCS) for (i, disp) in zip(inc, eachcol(u))]...))))

    undisp = Point3.([p.nodeStart.position, p.nodeEnd.position])
end
begin
    lw = 5
    fig = Figure(backgroundcolor = :black)
    ax = Axis3(fig[1,1], 
        aspect = :data,
        # aspect = (1,1,1)
        )

    # hidedecorations!(ax)
    simplifyspines!(ax)

    linesegments!(undisp,
    color = :white,
    linestyle = :dash)

    lines!(xyz,
        linewidth = lw,
        color = :white)

    on(df) do _ reset_limits!(ax) end

    fig
end

for i = 1:100
    df[] = i
    sleep(1e-3)
end

uglobal = [p.nodeStart.displacement; p.nodeEnd.displacement]
ulocal = p.R * uglobal
L = p.length

xrange = range(0, L, n)

uX = ulocal[[1, 7]]
uY = ulocal[[2,6,8,12]]
uZ = ulocal[[3,5,9,11]]

hcat([AsapToolkit.Naxial(x, L) * uX for x in xrange]...)
hcat([AsapToolkit.N(x, L) * uY for x in xrange]...)
hcat([AsapToolkit.N(x,L) * uZ for x in xrange]...)