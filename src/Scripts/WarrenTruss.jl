# regular warren truss

begin
    n = 11
    dx = 1500
    dy = 1400

    section = toASAPtruss(rand(allHSSRound()), Steel_Nmm.E)
end

warrentruss = generatewarren2d(n,
    dx,
    dy,
    section)

truss = warrentruss.model

begin
    pts = Point3.(getproperty.(truss.nodes, :position))
    els = vcat([pts[e.nodeIDs] for e in truss.elements]...)
    fs = getindex.(getproperty.(truss.elements, :forces), 2)
    cr = maximum(abs.(fs)) .* (-1, 1) .* .75
end

begin
    fig = Figure(resolution = (1000,500))
    ax = Axis(fig[1,1],
        aspect = DataAspect())

    hidedecorations!(ax); hidespines!(ax)

    linesegments!(els,
        linewidth = 5,
        color = fs,
        colormap = pink2blue,
        colorrange = cr)

    cb = Colorbar(fig[2,1],
        vertical = false,
        flipaxis = false,
        colorrange = cr ./ .75,
        colormap = pink2blue,
        label = "P [N]")

    fig
end


#irregular
xpos = cumsum(rand(7) .* 2000) .+ 500
ypos = - sort(rand(7) .* 1000); ypos[1] = ypos[end] = 0.
ypos2 = sort(rand(6) .* 1500) .+ 1500

truss = generatewarren2d(xpos,
    ypos,
    ypos2,
    section)

begin
    pts = Point3.(getproperty.(truss.nodes, :position))
    els = vcat([pts[e.nodeIDs] for e in truss.elements]...)
    fs = getindex.(getproperty.(truss.elements, :forces), 2)
    cr = maximum(abs.(fs)) .* (-1, 1) .* .75
end

begin
    fig = Figure(resolution = (1000,500))
    ax = Axis(fig[1,1],
        aspect = DataAspect())

    hidedecorations!(ax); hidespines!(ax)

    linesegments!(els,
        linewidth = 5,
        color = fs,
        colormap = pink2blue,
        colorrange = cr)

    cb = Colorbar(fig[2,1],
        vertical = false,
        flipaxis = false,
        colorrange = cr ./ .75,
        colormap = pink2blue,
        label = "P [N]")

    fig
end
