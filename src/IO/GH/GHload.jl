abstract type GHload end
function GHload end

struct GHnodeforce <: GHload
    iNode::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::NodeForce)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodeforce(i, value, id)
end

struct GHnodemoment <: GHload
    iNode::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::NodeMoment)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodemoment(i, value, id)

end

struct GHlineload <: GHload
    iElement::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::LineLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHlineload(i, value, id)

end

struct GHpointload <: GHload
    iElement::Int64
    x::Float64
    value::Vector{Float64}
    id::String
end

function GHload(load::PointLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)
    x = load.position

    return GHpointload(i, x, value, id)

end
