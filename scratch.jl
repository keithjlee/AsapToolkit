using Asap, AsapToolkit, kjlMakie
set_theme!(kjl_dark)
#sections
begin
    W = allW()
    iStiffness = sortperm(getproperty.(W, :Ix))
    nsecs = length(W)

    mat = Steel_Nmm
    girder = toASAPframe((W[iStiffness])[Int(round(nsecs * 0.5))], mat.E, mat.G)
    joist = toASAPframe((W[iStiffness])[Int(round(nsecs * .2))], mat.E, mat.G)

    iArea = sortperm(getproperty.(W, :A))
    col = toASAPframe((W[iArea])[Int(round(nsecs * .75))], mat.E, mat.G)

    H = allHSSRound()
    iArea = sortperm(getproperty.(H, :A))
    nh = length(H)
    tube = toASAPframe((H[iArea])[Int(round(nh * 0.9))], mat.E, mat.G)

    girder.ρ = joist.ρ = tube.ρ = mat.ρ
end

begin
    begin
        nx = 6
        ny = 6
        nz = 50
        dx = 4000
        dy = 3200
        dz = 3500
    end

    @time frame = AsapToolkit.generateframe(nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        1000,
        col,
        girder,
        joist,
        tube;
        # joistPsi = 0,
        columnPsi = pi/2,
        joistRelease = :joist,
        primaryRelease = :fixedfixed);

    model = frame.model;
    loads = [LineLoad(j, [0., 0., -10.]) for j in model.elements[:joist]];
    @time solve!(model, loads)

    dispfac = Observable(2.)
    resolution = 20

    ixn = frame.iExteriorXnodes
    iyn = frame.iExteriorYnodes
    ixj = frame.iExteriorXjoists
    iyp = frame.iExteriorYprimaries

    windX = [LineLoad(j, [10., 0., 0.]) for j in model.elements[ixj]]
    windY = [LineLoad(j, [0., 10., 0.]) for j in model.elements[iyp]]
    @time solve!(model, windY);

    moms = vcat([e.forces[[6,12]] for e in model.elements]...)
    momrange = maximum(abs.(moms)) .* (-1,1)

    begin
        N1 = Point3.([n.position for n in model.nodes])
        E1 = vcat([N1[e.nodeIDs] for e in model.elements]...)

        N2 = @lift([Point3(n.position .+ $dispfac * n.displacement[1:3]) for n in model.nodes])
        E2 = [@lift(Point3.(eachcol(displacedshape(e; factor = $dispfac, n= resolution)))) for e in model.elements]
        E2simple = @lift(vcat([$N2[e.nodeIDs] for e in model.elements]...))

        axf = getindex.(getproperty.(model.elements, :forces), 7)
        cr = maximum(abs.(axf)) .* (-1,1)
        lw = abs.(axf ./ maximum(abs.(axf))) .* 10  #.+ 1
        lw2 = abs.(moms ./ maximum(abs.(moms))) .* 10  #.+ 1

        nx = @lift($N2[ixn])
        ny = @lift($N2[iyn])
        jx = @lift(vcat([$N2[e.nodeIDs] for e in model.elements[ixj]]...))
        py = @lift(vcat([$N2[e.nodeIDs] for e in model.elements[iyp]]...))


        axf2 = @lift(axf .* $dispfac)
        moms2 = @lift(moms .* $dispfac)
        
        fig = Figure(backgroundcolor = :black, resolution = (1000,1000))
        ax = Axis3(fig[1,1],
            aspect = :data,
            # aspect = (1,1,1)
            )

        # simplifyspines!(ax)
        # gridtoggle!(ax)
        # hidespines!(ax)
        # gridtoggle!(ax)

        hidedecorations!(ax)
        hidespines!(ax)

        e_undisp = linesegments!(E1,
            color = :white,
            # linestyle = :dash,
            linewidth = .2)

        # scatter!(N1)
        
        # e_disp = lines!.(E2,
        #     color = :white,
        #     )

        e_simp = linesegments!(E2simple,
            # color = sign.(axf),
            color = axf2,
            colorrange = cr,
            colormap = pink2blue,
            linewidth = lw,
            )

        # e_simp_mom = linesegments!(E2simple,
        #     # color = sign.(axf),
        #     color = moms2,
        #     colorrange = momrange,
        #     colormap = pink2blue,
        #     linewidth = lw2,
        #     )

        pnx = scatter!(nx,
            color = green)
        pny = scatter!(ny,
            color = green)
        ljx = linesegments!(jx,
            color = green,
            linewidth = 5)
        lpy = linesegments!(py,
            color = green,
            linewidth = 5)

        pnx.visible = pny.visible = ljx.visible = lpy.visible = false

        
        on(dispfac) do _
            reset_limits!(ax)
        end

        fig
    end
end



begin
    n1 = Node([0., 0., 0.], :fixed)
    n2 = Node([5e3, 0., 0.], :fixed)

    nodes = [n1, n2]

    e = Element(n1, n2, joist, :joist)
    elements = [e]

    loads = [LineLoad(e, [0., 0., -20.])]

    m2 = Model(nodes, elements, loads)
    solve!(m2)
end



#single order
begin
    l = loads[1]

    xrange = range(0, e.length, 100)
    w = 20.

    mfunc = AsapToolkit.MfuncDict[(typeof(l), e.release)]
    vfunc = AsapToolkit.VfuncDict[(typeof(l), e.release)]
    dfunc = AsapToolkit.DfuncDict[(typeof(l), e.release)]

    m = mfunc.(w, e.length, xrange)
    v = vfunc.(w, e.length, xrange)
    d = dfunc.(w, e.length, xrange, e.section.E, e.section.Izz)
end
##### truss
begin
    n = 25
    dx = 3000
    dy = 6000
    section = toASAPframe(rand(allW()), mat.E, mat.G)
end

truss = generatewarren2d(n, dx, dy, section)

model = truss.model
L1 = [NodeForce(n, [0., -50e3, 0.]) for n in model.nodes[:bottomchord]]
@time solve!(model, L1);

begin
    pts = Point3.([n.position for n in model.nodes])
    els = vcat([pts[e.nodeIDs] for e in model.elements]...)
    axf = last.(getproperty.(model.elements, :forces))
    cr = maximum(abs.(axf)) .* (-1, 1)
    lw = abs.(axf) ./ maximum(abs.(axf)) .* 10

    dfac = Observable(2.)

    dpts = @lift([Point3(n.position .+ $dfac * n.displacement) for n in model.nodes])
    dels = @lift(vcat([$dpts[e.nodeIDs] for e in model.elements]...))

    fig = Figure(backgroundcolor = :black)
    ax = Axis(fig[1,1],
        aspect = DataAspect())

    hidespines!(ax)

    e2 = linesegments!(dels,
        color = axf,
        colorrange = cr,
        # linewidth = lw,
        linewidth = 10,
        colormap = pink2blue)

    on(dfac) do _
        reset_limits!(ax)
    end

    fig
end

for i = 0:.1:500
    dfac[] = i
    sleep(1e-3)
end