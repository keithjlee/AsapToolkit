begin
    
    # n1 = Node([0., 0., 0.], :pinned)
    # n2 = Node([7500., 1000., 0.], :pinned)

    # n3 = Node([-500., 400., 500.], n1.dof)
    # n4 = Node([0., 6500., 3000.], n2.dof)

    n1 = Node([0., 0., 0.], :fixed)
    n2 = Node([0., 7500., 0.], :fixed)

    n3 = Node([5000., 0., 0.], n1.dof)
    n4 = Node([5000., 7500., 0.], n2.dof)

    nodes = [n1, n2, n3, n4]

    e1 = Element(n1, n2, girder)
    e2 = Element(n3, n4, girder)

    joists = [BridgeElement(e1, x, e2, x, joist, :joist) for x in range(0, 1, 10)[2:end-1]]

    for el in joists el.id = :joist; end

    elements = [e1, e2, joists...]

    loads = [LineLoad(j, [0., 0., -10]) for j in joists]

    model = Model(nodes, elements, loads)
    solve!(model)

end


ploads = [NodeForce(n, [-5e3, 0., 45e3]) for n in model.nodes[5:12]]
ploads2 = [NodeForce(n, [2.5e3, 0., 0.]) for n in model.nodes[13:20]]
ploads3 = [PointLoad(j, rand(), [0., 0., -20e3]) for j in model.elements[:joist]]
combined = [ploads; model.loads; ploads2; ploads3]
solve!(model, combined)

# for j in model.elements[:joist] j.release = :fixedfixed end

begin
    dfac = Observable(100.)
    resolution = 20

    n_ud = Point3.(getproperty.(model.nodes, :position))
    e_ud = vcat([n_ud[e.nodeIDs] for e in model.elements]...)

    n_d = @lift(n_ud .+ $dfac * [n.displacement[1:3] for n in model.nodes])

    e_d = [@lift(displacedshape(e; factor = $dfac, n = resolution)) for e in model.elements]

    mids = Point3.(vcat([repeat([midpoint(e)], 3) for e in model.elements[:joist]]...))
    lcs = 500 .* Vec3.(vcat(getproperty.(model.elements[:joist], :LCS)...))
    cols = repeat([pink, blue, green], length(model.elements[:joist]))

    begin
        fig = Figure(backgroundcolor = :black, resolution = (2000,2000))
        ax = Axis3(fig[1,1],
            aspect = :data,
            zlabeloffset = 100,
            protrusions = 200)

        labelscale!(ax, 2)

        ud = linesegments!(e_ud,
            color = :white,
            linestyle = :dash)

        l = arrows!(mids, lcs,
            color = cols,
            linewidth = 50,
            arrowsize = 75)

        d = [lines!(ed,
            color = :white,
            linewidth = 3) for ed in e_d]

        tt = GLMakie.text!(n_ud, text = string.(getproperty.(model.nodes, :nodeID)),
            fontsize = 40)

        on(dfac) do _
            reset_limits!(ax)
        end

        fig
    end
end

j = rand(model.elements)

res = internalforces(j, model)

results = [internalforces(e, model) for e in model.elements]

#find shattered elements
ichecked = Vector{Int64}()
inds = Vector{Vector{Int64}}()
eids = getproperty.(model.elements, :elementID)
for e in model.elements
    id = e.elementID

    in(id, ichecked) && continue

    igroup = findall(eids .== id)

    push!(inds, igroup)
    push!(ichecked, id)
end

results = [internalforces(model.elements[id], model) for id in inds]

xvals = getproperty.(results, :x)
myvals = getproperty.(results, :My)
vyvals = getproperty.(results, :Vy)
mzvals = getproperty.(results, :Mz)
vzvals = getproperty.(results, :Vz)

begin
    fig = Figure(backgroundcolor = :black)
    axM = Axis(fig[1,1],
        yreversed = true,
        aspect = nothing)

    axV = Axis(fig[2,1],
        aspect = nothing)

    axMz = Axis(fig[1,2],
        yreversed = true,
        aspect = nothing)

    axVz = Axis(fig[2,2],
        aspect = nothing)

    hlines!.((axM, axV, axMz, axVz), [0.], color = :white)

    lines!.(axM, xvals, myvals,
        color = (blue, 0.5),
        linewidth = 4)

    lines!.(axV, xvals, vyvals,
        color = (green, 0.5),
        linewidth = 4)

    lines!.(axMz, xvals, mzvals,
        color = (blue, 0.5),
        linewidth = 4)

    lines!.(axVz, xvals, vzvals,
        color = (green, 0.5),
        linewidth = 4)

    fig
end


#beam information
element = rand(model.elements[:joist])
release = element.release
L = element.length

#discretization
xinc = collect(range(0, L, resolution))

#end node information
uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]

# end forces that are relevant to the given element/release condition
Flocal = (element.R * element.K * uglobal) .* release2DOF[release]

# shear/moment acting at the *starting* point of an element in LCS
Vystart, Mystart, Vzstart, Mzstart = Flocal[[2, 6, 3, 5]] .* [1, 1, 1, -1]

# initialize internal force vectors
My = Vystart .* xinc .- Mystart
Vy = zero(My) .+ Vystart
Mz = Vzstart .* xinc .- Mzstart
Vz = zero(Mz) .+ Vzstart

# accumulate loads
for load in model.loads[element.loadIDs]
    AsapToolkit.accumulate!(load,
        xinc,
        My,
        Vy,
        Mz,
        Vz)  
end

begin
    begin
        fig = Figure(backgroundcolor = :black)
        axM = Axis(fig[1,1],
            yreversed = true,
            aspect = nothing)
    
        axV = Axis(fig[2,1],
            aspect = nothing)
    
        hlines!.((axM, axV), [0.], color = :white)
    
        lines!(axM, xinc, Mz,
            color = blue,
            linewidth = 4)
    
        lines!(axV, xinc, Vz,
            color = green,
            linewidth = 4)
    
        fig
    end
end