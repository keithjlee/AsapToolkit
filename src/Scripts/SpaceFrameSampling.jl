using Interpolations, kjlMakie, Asap, AsapToolkit, JSON;
set_theme!(kjl_dark)

nxrange = 10:30
nyrange = 10:30
dxrange = 1000:250:2000
dyrange = 1000:250:2000
dzrange = 1000:250:3500
supps = [:corner, :x, :y, :xy]

begin
    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

@time begin
    nx = rand(nxrange)
    ny = rand(nyrange)
    dx = rand(dxrange)
    dy = rand(dyrange)
    dz = rand(dzrange)
    supp = rand(supps)

    sf = generatespaceframe(
        nx, 
        dx, 
        ny, 
        dy, 
        dz, 
        tube, 
        ; 
        load = [0., 0., -30e3],
        support = supp);

    truss = sf.truss;

    # dispfac = Observable(5.)
    # crfac = Observable(0.5)
    # lw = Observable(3.)

    # sizes = trusssizer.(truss.elements, σ)
    # AI = Observable([Point2(size) for size in sizes if last(size) != 0])
    # A = Observable([first(size) for size in sizes if last(size) == 0])

    # eforces = Observable(getindex.(getproperty.(truss.elements, :forces), 2))
    # cr = @lift($crfac .* (-1, 1) .* maximum(abs.($eforces)))

    # pts = Observable(Point3.(getproperty.(truss.nodes, :position)))
    # disps = Observable(getproperty.(truss.nodes, :displacement))
    # els = @lift(vcat([$pts[id] for id in getproperty.(truss.elements, :nodeIDs)]...))

    # p_d = @lift($pts .+ $dispfac .* $disps)
    # e_d = @lift(vcat([$p_d[id] for id in getproperty.(truss.elements, :nodeIDs)]...))

    # p_supports = Observable(Point3.(getproperty.(truss.nodes[:support], :position)))

    # fig = Figure(resolution = (1000,500))

    # axtop = Axis3(fig[1,1],
    #     aspect = :data,
    #     azimuth = pi,
    #     title = "Supports [mm]",
    #     titlesize = 15,
    #     zticksvisible = false,
    #     zticklabelsvisible = false,
    #     xlabelvisible = false,
    #     zlabelvisible = false,
    #     ylabelvisible = false,
    #     xticks = 0:1e4:2e4,
    #     yticks = 0:2e4:4e4,
    #     elevation = pi/2)

    # # hidedecorations!(axtop); hidespines!(axtop)
    # hidespines!(axtop); gridtoggle!(axtop)
    # axtop.titlevisible = true

    # xy = linesegments!(els,
    #     linewidth = lw,
    #     color = eforces,
    #     colorrange = cr,
    #     colormap = pink2blue,
    #     # color = :white
    #     )

    # suppnodes = scatter!(p_supports,
    #     color  = :black,
    #     strokecolor = :white,
    #     markersize = 10,
    #     overdraw = true)

    # ax = Axis3(fig[1,2],
    #     aspect = :data)

    # gridtoggle!(ax); simplifyspines!(ax)

    # u = linesegments!(els,
    #     color = :white)

    # u.visible = false

    # d = linesegments!(e_d,
    #     color = eforces,
    #     colorrange = cr,
    #     colormap = pink2blue,
    #     linewidth = lw
    #     )

    # hidedecorations!(ax); hidespines!(ax)

    # axComp = Axis(fig[2,1],
    #     xlabel = "min. A [mm²]",
    #     ylabel = "min. I [mm⁴]",
    #     title = "Compression",
    #     aspect = nothing)

    # scatter!(AI, color = pink)

    # axTens = Axis(fig[2,2],
    #     xlabel = "min. A [mm²]",
    #     title = "Tension",
    #     aspect = nothing)

    # vlines!(A, color = (blue, 0.25))

    # labelscale!.((axComp, axTens), 0.75)

    # on(dispfac) do _
    #     reset_limits!(ax)
    # end


    # fig
end;
