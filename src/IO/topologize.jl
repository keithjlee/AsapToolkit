"""
    topologize(model::Model; one_based = false, supplementary_data = nothing)

Export the topology of an Asap model as a JSON file.

# Arguments
- `model::Model` model to export (must be solved)

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
function topologize(model::Model; one_based = false, supplementary_data = nothing)

    isnothing(model.results) && error("Analyze model before export")

    xyz = node_positions(model)

    i_starts = [element.nodeStart.index for element in model.elements]
    i_ends = [element.nodeEnd.index for element in model.elements]

    i_support = [i for i = 1:length(model.nodes) if !all(model.nodes[i].fixity)]

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
        "axialForces" => [axial_force(model.results, element) for element in model.elements],
        "areas" => getproperty.(getproperty.(model.elements, :section), :A)
    )

    if !isnothing(supplementary_data) && typeof(supplementary_data) <: Dict
        out_dict = merge(out_dict, supplementary_data)
    end

    return out_dict

end
