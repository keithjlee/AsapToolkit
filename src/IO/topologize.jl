function topologize(model::Asap.AbstractModel; one_based = false)

    xyz = node_positions(model)

    indices = getproperty.(model.elements, :nodeIDs)

    i_starts = getindex.(indices, 1)
    i_ends = getindex.(indices, 2)

    i_support = [i for i = 1:model.nNodes if !all(model.nodes[i].dof)]

    if !one_based
        i_starts .-= 1
        i_ends .-= 1
        i_support .-= 1
    end

    out_dict = Dict(
        "x" => xyz[:, 1],
        "y" => xyz[:, 2],
        "z" => xyz[:, 3],
        "iStart" => i_starts,
        "iEnd" => i_ends,
        "iNodes" => i_support,
        "axialForces" => axial_force.(model.elements),
        "areas" => getproperty.(getproperty.(model.elements, :section), :A)
    )

    return out_dict

end