using SteelSections, kjlMakie
set_theme!(kjl_light)

W = allW()

webArea(w::SteelSections.W) = (w.d - 2w.tf) * w.tw

relfields = (:d, :A, :bf, :tf, :tw, :Ix, :Zx, :Sx, :Iy, :Zy, :Sy, :J)
l = length(relfields)

begin
    fig = Figure()

    for i = 1:l
        for j = 1:l
            ax = Axis(fig[i, j],
                aspect = 1)

            hidedecorations!(ax)
            hidespines!(ax)

            scatter!(getproperty.(W, relfields[i]), getproperty.(W, relfields[j]),
                color = getproperty.(W, :A),
                markersize = 5,
                colormap = :tempo)
        end
    end

    fig
end


fac = 25
d = 3

x1 = rand(d, 20) .* fac
x1[1:2,:] = sort(x1[1:2,:], dims = 2)
x2 = rand(d, 31) .* fac
x2[1:2,:] = sort(x2[1:2,:], dims = 2)

x1 = [1. 2. 3. 4. 5.; 0. 4. 3.9 2. 3.1; 2. 4. 5. 6. 10.]
x2 = [2. 3. 5. 7.; 1. 5. 2.5 3.0; 3. 5. 5.5 9.]


L(x) = sum(norm.(eachcol(diff(x, dims = 2))))

# total length of curves
l1 = L(x1)
l2 = L(x2)

# parameterization
t1 = range(0, 1, size(x1)[2])
t2 = range(0, 1, size(x2)[2])



p1 = ndInterpolate(x1)
p2 = ndInterpolate(x2)

nsamples = 100

trange1 = range(0, 1, nsamples)
trange2 = range(0, 1, nsamples)

freemap = [norm(p1(t1) .- p2(t2)) for t1 in trange1, t2 in trange2]


t1 = Observable(0.)
t2 = Observable(0.)
connector = @lift(Point3.([p1($t1), p2($t2)]))
point2 = @lift(Point2($t1, $t2))
point3 = @lift(Point3($t1, $t2, norm(p1($t1) .- p2($t2))))
pointset2 = @lift([$point2])
pointset3 = @lift([$point3])

cm = Reverse(:grays);
cm = :magma;
begin
    fig = Figure()

    ax0 = Axis3(fig[1,1], aspect = (1,1,1))

    l1 = scatterlines!(Point3.(eachcol(x1)),
        color = blue)

    l2 = scatterlines!(Point3.(eachcol(x2)),
        color = green)

    c = lines!(connector,
        color = :black)

    ax = Axis(fig[1,2], aspect = nothing)

    heatmap!(trange1, trange2, freemap, interpolate = true, colormap = cm)
    # contour!(trange1, trange2, freemap, colormap = Reverse(:grays), overdraw = true, linewidth = 3)
    # vlines!()

    scatter!(point2, color = :white, strokecolor = :black)

    ax2 = Axis3(fig[1,3], aspect = (1,1,1))
    simplifyspines!(ax2)
    gridtoggle!(ax2)

    surface!(trange1, trange2, freemap, colormap = cm)
    contour3d!(trange1, trange2, freemap, colormap = Reverse(:grays), overdraw = true, linewidth = 3)
    scatter!(point3, color = :white, overdraw = true, strokecolor = :black)

    colsize!(fig.layout, 1, Aspect(1,1.))
    colsize!(fig.layout, 2, Aspect(1,1.))
    colsize!(fig.layout, 3, Aspect(1,1.))

    resize_to_layout!(fig)

    fig
end

rset1 = sort(rand(1000))

r2 = sort(abs.(randn(1000)))
rset2 = r2 ./ maximum(r2)

for (i,j) in zip(range(0, 1, 1000), range(.3, .89, 1000))
    t1[] = i
    t2[] = j

    sleep(.005)
end

for i = range(0, 1, 1000)
    t1[] = i
    t2[] = i
    sleep(.005)
end

M_udl(w, x, l) = w * x / 2 * (l - x)
V_udl(w,x, l) = w * (l /2 - x)

function M_4pt(P, x, l, a)
    if x < a
        P * x
    elseif a < x < l - a
        P * a
    else
        P * a - (P * (x - l + a))
    end
end

function V_4pt(P, x, l, a)
    if x < a
        P
    elseif a < x < l - a
        0.
    else
        -P
    end
end

l = 12
a = 3
P = 34
range1 = range(0, l, 200)
x1 = Float64.([(collect(range1) .- l/2)'; M_4pt.(P, range1, l, a)'; V_4pt.(P, range1, l, a)'])
# x1[3,:] .= abs.(x1[3,:])

l2 = 8
w = 25
range2 = range(0, l2, 200)
x2 = Float64.([(collect(range2) .- l2/2)'; M_udl.(w, range2, l2)'; V_udl.(w, range2, l2)'])
# x2[3,:] .= abs.(x2[3,:])

d1 = collect(range1) .- l/2
d2 = collect(range2) .- l2/2

p1 = ndInterpolate(x1, d1)
p2 = ndInterpolate(x2, d2)

nsamples = 100

trange1 = d1
trange2 = d2

freemap = [norm(p1(t1) .- p2(t2)) for t1 in trange1, t2 in trange2]


t1 = Observable(0.)
t2 = Observable(0.)
connector = @lift(Point3.([p1($t1), p2($t2)]))
point2 = @lift(Point2($t1, $t2))
point3 = @lift(Point3($t1, $t2, norm(p1($t1) .- p2($t2))))
pointset2 = @lift([$point2])
pointset3 = @lift([$point3])

xloss = Observable(Vector{Float64}())
loss = Observable(Vector{Float64}())

cm = Reverse(:grays);
cm = :magma; set_theme!(kjl_dark)
cm = :tempo;
begin
    fig = Figure(backgroundcolor = :black)

    ax0 = Axis3(fig[1,1], aspect = (1,1,1), xlabel = "x [m]",
        ylabel = "M [kNm]",
        zlabel = "V [kN]")

    ax0.protrusions = 50

    simplifyspines!(ax0)
    gridtoggle!(ax0)

    l_1 = scatterlines!(Point3.(eachcol(x1)),
        markersize = 2,
        linewidth = 3,
        color = blue)

    l_2 = scatterlines!(Point3.(eachcol(x2)),
        markersize = 2,
        linewidth = 3,
        color = green)

    c = scatterlines!(connector,
        color = :white)

    ax = Axis(fig[1,2], aspect = nothing,
        xlabel = "4 pt load",
        xlabelcolor = blue,
        xtickcolor = blue,
        xticklabelcolor = blue,
        ylabel = "UDL",
        ylabelcolor = green,
        ytickcolor = green,
        yticklabelcolor = green
        )

    heatmap!(trange1, trange2, freemap, interpolate = true, colormap = cm)
    contour!(trange1, trange2, freemap, colormap = Reverse(:grays), overdraw = true, linewidth = 3, levels = 10)
    # vlines!()

    scatter!(point2, color = :white, strokecolor = :black)

    ax2 = Axis3(fig[1,3], aspect = (1,1,1),
        xlabel = "4 pt load [m]",
        xlabelcolor = blue,
        xtickcolor = blue,
        xticklabelcolor = blue,
        ylabel = "UDL [m]",
        ylabelcolor = green,
        ytickcolor = green,
        yticklabelcolor = green,
        zlabel = "||VM(xUDL) - VM(x4pt)||")
    ax2.protrusions = 50
    simplifyspines!(ax2)
    gridtoggle!(ax2)

    surface!(trange1, trange2, freemap, colormap = cm)
    # contour3d!(trange1, trange2, freemap, colormap = Reverse(:grays), overdraw = true, linewidth = 3)
    scatter!(point3, color = :white, overdraw = true, strokecolor = :black)

    axloss = Axis(fig[2,:],
        aspect = nothing,
        xlabel = "x [m]",
        ylabel = "Distance")

    lines!(xloss, loss, color = :white)

    on(t1) do _
        reset_limits!(axloss)
    end

    colsize!(fig.layout, 1, Aspect(1,1.))
    colsize!(fig.layout, 2, Aspect(1,1.))
    colsize!(fig.layout, 3, Aspect(1,1.))

    resize_to_layout!(fig)

    fig
end

xloss = Observable(Vector{Float64}())
loss = Observable(Vector{Float64}())

for i in d2
    t1[] = i
    t2[] = i

    push!(xloss[], t2[])
    push!(loss[], point3[][3])

    notify.((xloss, loss))
    sleep(.01)
end

iterator = range(0, 1, 600)
anginc = pi / length(iterator)

iterator = range(-4, 4, length = 300)
anginc = pi / length(iterator)
record(fig, "frechet_dark.gif", iterator; framerate = 30) do x
    ax0.azimuth[] += anginc
    ax2.azimuth[] += anginc
    t1[] = x
    t2[] = x

    push!(xloss[], t2[])
    push!(loss[], point3[][3])

    notify.((xloss, loss))
end