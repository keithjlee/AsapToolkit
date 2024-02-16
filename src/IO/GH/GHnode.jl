struct GHnode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    reaction::Vector{Float64}
    u::Vector{Float64}
    displacement::Vector{Float64}
    id::String

    function GHnode(node::TrussNode)
        position = node.position
        dof = [node.dof; [true, true, true]]
        nodeID = node.nodeID - 1
        reaction = node.reaction
        u = node.displacement
        displacement = node.displacement[1:3]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end

    function GHnode(node::Node)
        position = node.position
        dof = node.dof
        nodeID = node.nodeID - 1
        reaction = node.reaction
        u = node.displacement
        displacement = node.displacement[1:3]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end
end