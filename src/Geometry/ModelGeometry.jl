struct ModelGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    forces::Vector{Float64}
    max_abs_force::Float64
    moments::Vector{Vector{Float64}}
    max_abs_moment::Float64
    areas::Vector{Float64}
    max_area::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}

    function ModelGeo(model::Model)

        nodes = getproperty.(model.nodes, :position)
        nodes_xy = [node[1:2] for node in nodes]

        disp = getproperty.(model.nodes, :displacement)
        disp_xy = [d[1:2] for d in disp]

        indices = getproperty.(model.elements, :nodeIDs)
        forces = getindex.(getproperty.(model.elements, :forces), 7)
        max_abs_force = maximum(abs.(forces))

        moments = [element.forces[[4,5,6,10,11,12]] for element in model.elements]
        moment_magnitudes = vcat([[norm(moment[1:3]), norm(moment[4:end])] for moment in moments]...)
        max_abs_moment = maximum(moment_magnitudes)

        areas = getproperty.(getproperty.(model.elements, :section), :A)
        max_area = maximum(areas)

        element_vectors = Asap.localx.(model.elements)
        element_vectors_xy = [evec[1:2] for evec in element_vectors]

        return new(
            nodes,
            nodes_xy,
            disp,
            disp_xy,
            indices,
            forces,
            max_abs_force,
            moments,
            max_abs_moment,
            areas,
            max_area,
            getproperty.(model.elements, :length),
            element_vectors,
            element_vectors_xy
        )
    end

end