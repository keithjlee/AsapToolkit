using Asap, AsapToolkit, kjlMakie
set_theme!(kjl_dark)
#sections
tube = toASAPtruss(rand(allHSSRound()), Steel_Nmm.E)

begin
    nx = 20
    dx = 1000.
    ny = 31
    dy = 1300.
    dz = 1250.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

@time sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3]);
truss = sf.truss;

#randomly fix nodes



#plot objects
begin
    dispfac = Observable(5.)
    crfac = Observable(0.5)

    sizes = trusssizer.(truss.elements, σ)
    AI = Observable([Point2(size) for size in sizes if last(size) != 0])
    A = Observable([first(size) for size in sizes if last(size) == 0])

    eforces = Observable(getindex.(getproperty.(truss.elements, :forces), 2))
    cr = @lift($crfac .* (-1, 1) .* maximum(abs.($eforces)))

    pts = Observable(Point3.(getproperty.(truss.nodes, :position)))
    disps = Observable(getproperty.(truss.nodes, :displacement))
    els = @lift(vcat([$pts[id] for id in getproperty.(truss.elements, :nodeIDs)]...))

    p_d = @lift($pts .+ $dispfac .* $disps)
    e_d = @lift(vcat([$p_d[id] for id in getproperty.(truss.elements, :nodeIDs)]...))

    p_supports = Observable(Point3.(getproperty.(truss.nodes[:support], :position)))
end

begin
    fig = Figure(resolution = (1000,500))

    axtop = Axis3(fig[1,1],
        aspect = :data,
        azimuth = pi,
        title = "Supports [mm]",
        titlesize = 15,
        zticksvisible = false,
        zticklabelsvisible = false,
        xlabelvisible = false,
        zlabelvisible = false,
        ylabelvisible = false,
        xticks = 0:1e4:2e4,
        yticks = 0:2e4:4e4,
        elevation = pi/2)

    # hidedecorations!(axtop); hidespines!(axtop)
    hidespines!(axtop); gridtoggle!(axtop)
    axtop.titlevisible = true

    xy = linesegments!(els,
        linewidth = 2,
        color = eforces,
        colorrange = cr,
        colormap = pink2blue,
        # color = :white
        )

    supps = scatter!(p_supports,
        color  = :black,
        strokecolor = :white,
        markersize = 10,
        overdraw = true)

    ax = Axis3(fig[1,2],
        aspect = :data)

    gridtoggle!(ax); simplifyspines!(ax)

    u = linesegments!(els,
        color = :white)

    u.visible = false

    d = linesegments!(e_d,
        color = eforces,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 2
        )

    hidedecorations!(ax); hidespines!(ax)

    axComp = Axis(fig[2,1],
        xlabel = "min. A [mm²]",
        ylabel = "min. I [mm⁴]",
        title = "Compression",
        aspect = nothing)

    scatter!(AI, color = pink)

    axTens = Axis(fig[2,2],
        xlabel = "min. A [mm²]",
        title = "Tension",
        aspect = nothing)

    vlines!(A, color = (blue, 0.25))

    labelscale!.((axComp, axTens), 0.75)

    on(dispfac) do _
        reset_limits!(ax)
    end


    fig
end

#perimeter squares
begin
    x_perimeter1 = sf.isquares[1,:]
    x_perimeter2 = sf.isquares[end,:]
    y_perimeter1 = sf.isquares[:,1]
    y_perimeter2 = sf.isquares[:,end]
end
nsquares = 2

# iterator = 1:400
# record(fig, "spaceframes.mp4", iterator; framerate = 20) do _
begin
    fixnode!.(truss.nodes, :free)
    for node in truss.nodes[:support]
        node.id = :free
    end

    x1set = rand(x_perimeter1, nsquares)
    for i in vcat(x1set...)
        fixnode!(truss.nodes[i], :pinned)
        truss.nodes[i].id = :support
    end

    x2set = rand(x_perimeter2, nsquares)
    for i in vcat(x2set...)
        fixnode!(truss.nodes[i], :pinned)
        truss.nodes[i].id = :support
    end

    y1set = rand(y_perimeter1, nsquares)
    for i in vcat(y1set...)
        fixnode!(truss.nodes[i], :pinned)
        truss.nodes[i].id = :support
    end

    y2set = rand(y_perimeter2, nsquares)
    for i in vcat(y2set...)
        fixnode!(truss.nodes[i], :pinned)
        truss.nodes[i].id = :support
    end

    updateDOF!(truss); solve!(truss)

    pts[] = Point3.(getproperty.(truss.nodes, :position))
    disps[] = getproperty.(truss.nodes, :displacement)
    p_supports[] = Point3.(getproperty.(truss.nodes[:support], :position))
    eforces[] = getindex.(getproperty.(truss.elements, :forces), 2)
    sizes = trusssizer.(truss.elements, σ)
    AI[] = [Point2(size) for size in sizes if last(size) != 0]
    A[] = [first(size) for size in sizes if last(size) == 0]

    reset_limits!.((ax, axComp, axTens))
end