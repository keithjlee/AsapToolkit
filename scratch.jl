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

    H = allHSSRound()
    iArea = sortperm(getproperty.(H, :A))
    nh = length(H)
    tube = toASAPframe((H[iArea])[Int(round(nh * 0.75))], mat.E, mat.G)

    girder.ρ = joist.ρ = tube.ρ = mat.ρ
end

nx = 4
ny = 3
nz = 10
dx = 4000
dy = 3200
dz = 3500

model = AsapToolkit.generateFrame(nx,
    dx,
    ny,
    dy,
    nz,
    dz,
    1000,
    girder,
    girder,
    joist,
    tube;
    joistRelease = :fixedfixed,
    primaryRelease = :freefree);

loads = [LineLoad(j, [0., 0., -20.]) for j in model.elements[:joist]]

solve!(model, loads)

dispfac = 200
resolution = 20
begin
    N1 = Point3.([n.position for n in model.nodes])
    E1 = vcat([N1[e.nodeIDs] for e in model.elements]...)

    N2 = [Point3(n.position .+ dispfac * n.displacement[1:3]) for n in model.nodes]
    E2 = [Point3.(eachcol(displacedshape(e; factor = dispfac, n= resolution))) for e in model.elements]
    E2simple = vcat([N2[e.nodeIDs] for e in model.elements]...)

    axf = getindex.(getproperty.(model.elements, :forces), 7)
    cr = maximum(abs.(axf)) .* (-1,1)
    lw = abs.(axf ./ maximum(abs.(axf))) .* 10
    
    
    fig = Figure(backgroundcolor = :black)
    ax = Axis3(fig[1,1],
        aspect = :data,
        # aspect = (1,1,1)
        )

    simplifyspines!(ax)
    gridtoggle!(ax)

    e_undisp = linesegments!(E1,
        color = :white,
        linestyle = :dash,
        linewidth = .5)

    # scatter!(N1)
    
    e_disp = lines!.(E2,
        color = :white,
        )

    # e_simp = linesegments!(E2simple,
    #     color = axf,
    #     colorrange = cr,
    #     colormap = pink2blue,
    #     # linewidth = lw,
    #     linewidth = 5,
    #     )


    fig
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

    disp = Point3.(eachcol(displacedshape(e; factor = 100)))
    lines(disp, axis = (type = Axis3, aspect = (1,1,1),))
end

begin
    E = mat.E
    I = joist.Izz

    L = 6000.
    n = 100
    x = range(0, L, n)

    P = 100e3 #N
    frac = 0.6
    w = 40. #N/mm=kN/m
end

#test Line load
begin
    m = MLine_freefree.(w, L, x)
    v = VLine_freefree.(w, L, x)
    d = DLine_freefree.(w, L, x, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

begin
    m = MLine_fixedfree.(w, L, x)
    v = VLine_fixedfree.(w, L, x)
    d = DLine_fixedfree.(w, L, x, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    hlines!(axV, [0.], color = :white)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    hlines!(axV, [0.], color = :white)
    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

begin
    m = MLine_freefixed.(w, L, x)
    v = VLine_freefixed.(w, L, x)
    d = DLine_freefixed.(w, L, x, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    hlines!(axV, [0.], color = :white)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    hlines!(axV, [0.], color = :white)
    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

begin
    m = MLine_fixedfixed.(w, L, x)
    v = VLine_fixedfixed.(w, L, x)
    d = DLine_fixedfixed.(w, L, x, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    hlines!(axV, [0.], color = :white)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    hlines!(axV, [0.], color = :white)
    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

##### Point Load
begin
    m = MPoint_freefree.(P, L, x, frac)
    v = VPoint_freefree.(P, L, x, frac)
    d = DPoint_freefree.(P, L, x, frac, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

begin
    m = MPoint_fixedfixed.(P, L, x, frac)
    v = VPoint_fixedfixed.(P, L, x, frac)
    d = DPoint_fixedfixed.(P, L, x, frac, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end

##combined
begin
    m = MPoint_freefree.(P, L, x, frac) .+ MLine_freefree.(w, L, x) .+ MPoint_freefree.(50e3, L, x, .25)
    v = VPoint_freefree.(P, L, x, frac) .+ VLine_freefree.(w, L, x) .+ VPoint_freefree.(50e3, L, x, .25)
    d = DPoint_freefree.(P, L, x, frac, E, I) .+ DLine_freefree.(w, L, x, E, I) .+ DPoint_freefree.(50e3, L, x, .25, E, I)

    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        aspect = nothing,
        yreversed = true,
        title = "M")

    hidexdecorations!(axM)
    lines!(x, m)

    axV = Axis(fig[2,1],
        aspect = nothing,
        title = "V")

    hidexdecorations!(axV)

    hlines!(axV, [0.], color = :white)
    lines!(x, v)

    axD = Axis(fig[3,1],
        aspect = nothing,
        yreversed = true,
        title = "Δ")

    lines!(x, d)

    [rowsize!(fig.layout, i, Aspect(1, 0.3)) for i = 1:3]
    resize_to_layout!(fig)

    fig
end