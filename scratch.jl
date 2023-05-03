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
    nx = 6
    ny = 6
    nz = 50
    dx = 4000
    dy = 3200
    dz = 3500
end

@time frame = AsapToolkit.generateFrame(nx,
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
    # columnPsi = pi/2,
    joistRelease = :fixedfixed,
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

wave = sin.(0:.05:4pi)
ainc = pi / length(wave)

record(fig, "oscillate.mp4", wave; framerate = 30) do x
    dispfac[] = x * 20
    ax.azimuth[] += ainc
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

joists = model.elements[:joist]
jforces = getindex.(getproperty.(joists, :forces), 7)
columns = model.elements[:column]
cforces = getindex.(getproperty.(columns, :forces), 7)

begin
    fig = Figure(backgroundcolor = :black)
    ax = Axis(fig[1,1],
        xlabel = "Axial Force [kN]",
        aspect = 1)

    # J = hist!(jforces ./ 1e3,
    #     color = blue)

    C = hist!(cforces ./ 1e3,
        color = green)

    fig

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