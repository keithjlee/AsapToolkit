using Asap, AsapToolkit, kjlMakie, Interpolations;
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
    fig = Figure(resolution = (1000,1000))
    ax = Axis3(fig[1,1],
        aspect = :data)

    d = linesegments!(e_d,
        color = eforces,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 2
        )

    hidedecorations!(ax); hidespines!(ax)

    cb = Colorbar(fig[2,1],
        vertical = false,
        colorrange = cr,
        colormap = pink2blue,
        flipaxis = false,
        tellheight = true,
        label = "P [N]")

    fig
end


## generate curvy object
n = 4
x = range(0,1,n)
y = range(0,1,n)
z = rand(n,n) .* 4500

itp = cubic_spline_interpolation((x,y), z)

@time sf = generatespaceframe(nx, dx, ny, dy, dz, itp, tube, true; load = [0., 0., -30e3]);
truss = sf.truss;

#plot objects
begin
    dispfac = Observable(1.)
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
    fig = Figure(resolution = (1000,1000))
    ax = Axis3(fig[1,1],
        aspect = :data)

    d = linesegments!(e_d,
        color = eforces,
        colorrange = cr,
        colormap = pink2blue,
        linewidth = 2
        )

    hidedecorations!(ax); hidespines!(ax)

    cb = Colorbar(fig[2,1],
        vertical = false,
        colorrange = cr,
        colormap = pink2blue,
        flipaxis = false,
        tellheight = true,
        label = "P [N]")

    fig
end