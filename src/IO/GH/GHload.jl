abstract type GHload end
function GHload end

struct GHnodeforce <: GHload
    value::Vector{Float64}
    id::String
    iNode::Int64
end

function GHload(load::NodeForce)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodeforce(value, id, i)
end

struct GHnodemoment <: GHload
    value::Vector{Float64}
    id::String
    iNode::Int64
end

function GHload(load::NodeMoment)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHnodemoment(value, id, i)

end

struct GHlineload <: GHload
    value::Vector{Float64}
    id::String
    iElement::Int64
end

function GHload(load::LineLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)

    return GHlineload(value, id, i)

end

struct GHpointload <: GHload
    value::Vector{Float64}
    id::String
    iElement::Int64
    x::Float64
end

function GHload(load::PointLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : string(load.id)
    x = load.position

    return GHpointload(value, id, i, x)

end

function categorize_loads(loads::Vector{Asap.AbstractLoad})

end