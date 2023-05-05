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
        nx = 8
        ny = 16
        nz = 2
        dx = 3200
        dy = 5000
        dz = 3500
    end

    @time frame = AsapToolkit.generateframe(nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        500,
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
end

begin
    windX = [LineLoad(j, [10., 0., 0.]) for j in model.elements[ixj]]
    windY = [LineLoad(j, [0., 10., 0.]) for j in model.elements[iyp]]
    @time solve!(model, windY);

    moms = vcat([e.forces[[6,12]] for e in model.elements]...)
    momrange = maximum(abs.(moms)) .* (-1,1)


    N1 = Point3.([n.position for n in model.nodes])
    E1 = vcat([N1[e.nodeIDs] for e in model.elements]...)

    N2 = @lift([Point3(n.position .+ $dispfac * n.displacement[1:3]) for n in model.nodes])
    E2 = [@lift(Point3.(eachcol(displacedshape(e; factor = $dispfac, n= resolution)))) for e in model.elements]
    # C2Ax = [e.forces[7] for e in model.elements]
    # C2r = maximum(abs.(C2Ax)) .* (-1, 1)

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
    
    
end

begin
    fig = Figure(backgroundcolor = :black, resolution = (1000,1000))
    ax = Axis3(fig[1,1],
        aspect = :data,
        )

    # hidedecorations!(ax)
    # hidespines!(ax)

    e_undisp = linesegments!(E1,
        color = :white,
        linewidth = .2)

    e_simp = linesegments!(E2simple,
        color = axf2,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = lw,
        )

    # E = lines!.(E2,
    #     # color = sign.(C2Ax),
    #     # colormap = pink2blue,
    #     # colorrange = C2r,
    #     # linewidth = 5
    #     color = :white
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

## rapid sampling
nxrange = collect(3:8)
nyrange = collect(3:10)
nzrange = collect(1:30)
dxrange = collect(2000:500:6000)
dyrange = collect(3000:500:8000)
dzrange = collect(3500:500:5000)
joistspacerange = collect(500:50:1000)
# buffrange = collect(10e3:5e3:20e3)

e1store = Vector{Vector{Point3{Float64}}}()
e2store = Vector{Vector{Point3{Float64}}}()

axfstore = Vector{Vector{Float64}}()
crstore = Vector{Tuple{Float64, Float64}}()
lwstore = Vector{Vector{Float64}}()
modelstore = Vector{Asap.Model}()

begin
    xbase = 0.
    ybase = 0.
    ylargest = 0.

    num_x = 10
    num_y = 10
    dispfac = 2
    resolution = 10
    buffer = 15e3
end

for i = 1:num_x * num_y
    println("Iteration $i")
    nx = rand(nxrange)
    dx = rand(dxrange)
    ny = rand(nyrange)
    dy = rand(dyrange)
    nz = rand(nzrange)
    dz = rand(dzrange)
    js = rand(joistspacerange)
    # buffer = rand(buffrange)


    @time frame = AsapToolkit.generateframe(nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        js,
        col,
        girder,
        joist,
        tube;
        base = [xbase, ybase, 0.],
        # joistPsi = 0,
        columnPsi = pi/2,
        joistRelease = :joist,
        primaryRelease = :fixedfixed);

    model = frame.model;

    ixn = frame.iExteriorXnodes
    iyn = frame.iExteriorYnodes
    ixj = frame.iExteriorXjoists
    iyp = frame.iExteriorYprimaries

    # windX = [LineLoad(j, [10., 0., 0.]) for j in model.elements[ixj]]
    # windY = [LineLoad(j, [0., 10., 0.]) for j in model.elements[iyp]]
    loads = [LineLoad(j, [0., 0., -20.]) for j in model.elements[:joist]]
    loads2 = [NodeForce(n, [10e3, 0., 0.]) for n in model.nodes[iyn]]

    @time solve!(model, [loads; loads2]);

    push!(modelstore, model)

    N1 = Point3.([n.position for n in model.nodes])
    E1 = vcat([N1[e.nodeIDs] for e in model.elements]...)

    N2 = [Point3(n.position .+ dispfac * n.displacement[1:3]) for n in model.nodes]
    E2simple = vcat([N2[e.nodeIDs] for e in model.elements]...)

    axf = getindex.(getproperty.(model.elements, :forces), 7)
    cr = maximum(abs.(axf)) .* (-1,1)
    lw = abs.(axf ./ maximum(abs.(axf)))  #.+ 1


    push!(e1store, E1)
    push!(e2store, E2simple)
    push!(axfstore, axf)
    push!(crstore, cr)
    push!(lwstore, lw)

    xbase += nx * dx + buffer

    if ny * dy > ylargest
        ylargest = ny * dy
    end

    if i % num_x == 0
        ybase += ylargest + buffer
        xbase = 0
        ylargest = 0
    end

end

using CairoMakie; CairoMakie.activate!()
begin
    fig = Figure(backgroundcolor = :black, resolution = (6000,3000))

    ax = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax)
    hidespines!(ax)

    eu = [linesegments!(E, color = :white, linewidth = .5) for E in e1store]

    # ed = [linesegments!(E, 
    #     color = axf, 
    #     colormap = pink2blue, 
    #     colorrange = cr, 
    #     linewidth = 1
    #     ) for (E,axf,cr) in zip(e2store, axfstore, crstore)]

    fig
end

i = Observable(1)

E1 = @lift(e1store[$i])
E2 = @lift(e2store[$i])
axf = @lift(axfstore[$i])
cr = @lift(crstore[$i])
lw = @lift(4 * abs.($axf) ./ maximum(abs.($axf)) .+ 1)

anginc = 2pi / length(modelstore)
begin
    fig = Figure(backgroundcolor = :black, resolution = (1000,1000))

    ax1 = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax1)
    hidespines!(ax1)

    eu = linesegments!(E1, color = :white, linewidth = .5)


    ax2 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax2)
    hidespines!(ax2)

    eu = linesegments!(E2, 
        color = axf,
        colormap = pink2blue,
        colorrange = cr, 
        linewidth = lw)

    on(i) do _
        reset_limits!(ax1)
        reset_limits!(ax2)


    end

    fig
end

iterator = 1:100
record(fig, "highrises.mp4", iterator; framerate = 10) do x
    i[] = x
end