struct GHmodel
    nodes::Vector{GHnode}
    elements::Vector{GHelement}
    loads::Vector{GHload}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    dx::Vector{Float64}
    dy::Vector{Float64}
    dz::Vector{Float64}
    istart::Vector{Int64}
    iend::Vector{Int64}
    i_free_nodes::Vector{Int64}
    i_fixed_nodes::Vector{Int64}
end

"""
    GHmodel(model::Asap.AbstractModel)::GHmodel

Convert an Asap model into a condensed geometric format for interpoperability with other software.
"""
function GHmodel(model::Asap.AbstractModel)

    model.processed || (Asap.process!(model))

    nodes = GHnode.(model.nodes)
    elements = GHelement.(model.elements)
    loads = GHload.(model.loads)

    xyz = node_positions(model)

    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]

    istart = getproperty.(elements, :iStart)
    iend = getproperty.(elements, :iEnd)

    ifree = Vector{Int64}()
    ifixed = Vector{Int64}()

    dx = zero(x)
    dy = zero(y)
    dz = zero(z)

    for i in eachindex(nodes)

        if all(nodes[i].dof)
            push!(ifree, nodes[i].nodeID)
        else
            push!(ifixed, nodes[i].nodeID)
        end

        disp = nodes[i].displacement
        dx[i] = disp[1]
        dy[i] = disp[2]
        dz[i] = disp[3]
    end

    return GHmodel(
        nodes,
        elements,
        loads,
        x,
        y,
        z,
        dx,
        dy,
        dz,
        istart,
        iend,
        ifree,
        ifixed
    )

end
