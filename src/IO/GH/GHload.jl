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

const loadtype2vectorindex = Dict(
    GHnodeforce => 1,
    GHnodemoment => 2,
    GHlineload => 3,
    GHpointload => 4
)

function categorize_loads(loads::Vector{<:GHload})
    nodeforces = Vector{GHnodeforce}()
    nodemoments = Vector{GHnodemoment}()
    lineloads = Vector{GHlineload}()
    pointloads = Vector{GHpointload}()

    load_collectors = [nodeforces, nodemoments, lineloads, pointloads]

    for load in loads
        i = loadtype2vectorindex[typeof(load)]
        push!(load_collectors[i], load)
    end

    return load_collectors
end