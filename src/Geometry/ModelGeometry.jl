struct ModelGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    forces::Vector{Float64}
    max_abs_force::Float64
    Tx::Vector{Float64}
    max_abs_Tx::Float64
    My::Vector{Float64}
    max_abs_My::Float64
    Mz::Vector{Float64}
    max_abs_Mz::Float64
    areas::Vector{Float64}
    max_area::Float64
    Ix::Vector{Float64}
    max_Ix::Float64
    Iy::Vector{Float64}
    max_Iy::Float64
    J::Vector{Float64}
    max_J::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}

    function ModelGeo(model::Model)

        nodes = getproperty.(model.nodes, :position)
        nodes_xy = [node[1:2] for node in nodes]

        disp = [node.displacement[1:3] for node in model.nodes]
        disp_xy = [d[1:2] for d in disp]

        indices = getproperty.(model.elements, :nodeIDs)
        forces = getindex.(getproperty.(model.elements, :forces), 7)
        max_abs_force = maximum(abs.(forces))

        element_forces = getproperty.(model.elements, :forces)

        Tx = getindex.(element_forces, 10)
        max_abs_Tx = maximum(abs.(Tx))

        My = getindex.(element_forces, 11)
        max_abs_My = maximum(abs.(My))

        Mz = getindex.(element_forces, 12)
        max_abs_Mz = maximum(abs.(Mz))

        sections = getproperty.(model.elements, :section)

        areas = getproperty.(sections, :A)
        max_area = maximum(areas)

        Ix = getproperty.(sections, :Ix)
        max_Ix = maximum(Ix)

        Iy = getproperty.(sections, :Iy)
        max_Iy = maximum(Iy)

        J = getproperty.(sections, :J)
        max_J = maximum(J)

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
            Tx,
            max_abs_Tx,
            My,
            max_abs_My,
            Mz,
            max_abs_Mz,
            areas,
            max_area,
            Ix,
            max_Ix,
            Iy,
            max_Iy,
            J,
            max_J,
            getproperty.(model.elements, :length),
            element_vectors,
            element_vectors_xy
        )
    end

end