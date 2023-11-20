struct ModelGeo
    nodes::Vector{Vector{Float64}}
    nodes_xy::Vector{Vector{Float64}}
    disp::Vector{Vector{Float64}}
    disp_xy::Vector{Vector{Float64}}
    indices::Vector{Vector{Int64}}
    forces::Vector{Float64}
    max_abs_force::Float64
    moments::Vector{Float64}
    max_abs_moment::Float64
    lengths::Vector{Float64}
    element_vectors::Vector{Vector{Float64}}
    element_vectors_xy::Vector{Vector{Float64}}
    load_positions::Vector{Vector{Float64}}
    load_positions_xy::Vector{Vector{Float64}}
    load_vectors::Vector{Vector{Float64}}
    load_vectors_xy::Vector{Vector{Float64}}


end