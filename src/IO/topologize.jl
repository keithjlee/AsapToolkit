function topologize(model::Asap.AbstractModel; one_based = false)

    xyz = nodePositions(model)

    indices = getproperty.(model.elements, :nodeIDs)

    i_starts = getindex.(indices, 1)
    i_ends = getindex.(indices, 2)

    if one_based
        i_starts .+= 1
        i_ends .+= 1
    end

    i_support = [i for i = 1:model.nNodes if !all(model.nodes[i].dof)]

    out_dict = Dict(
        "x" => xyz[:, 1],
        "y" => xyz[:, 2],
        "z" => xyz[:, 3],
        "iStart" => i_starts,
        "iEnd" => i_ends,
        "iNodes" => i_support
    )

    return out_dict

end