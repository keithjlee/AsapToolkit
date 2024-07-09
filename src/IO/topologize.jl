"""
    topology(model::Asap.AbstractModel; one_based = false, supplementary_data = nothing)

Export the topology of an Asap model as a JSON file.

# Arguments
- `model::Asap.AbstractModel` model to export

## Optional Arguments
- `one_based::Bool = false` use 1-based indexing for topology
- `supplementary_data = nothing` a Dict object with additional data to export alongside the topology. E.g. if sections of elements are known, then one can do:
```julia
section_info = Dict(
    "widths" => section_widths,
    "depths" => section_depths
    )

topologize(model; supplementary_data = section_info)
```
"""
function topologize(model::Asap.AbstractModel; one_based = false, supplementary_data = nothing)

    xyz = node_positions(model)

    indices = Asap.nodeids(model.elements)

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

    if !isnothing(supplementary_data) && typeof(supplementary_data) <: Dict
        out_dict = merge(out_dict, supplementary_data)
    end

    return out_dict

end