abstract type AbstractParams end

struct TrussOptParams <: AbstractParams
    nodeids::Vector{Vector{Int64}} # [[iNodeStart, iNodeEnd] for element in elements]
    dofids::Vector{Vector{Int64}} # [[dofStartNode..., dofEndNode...] for element in elements]
    P::Vector{Float64} # External load vector
    freeids::Vector{Int64} # [DofFree1, DofFree2,...]
    inzs::Vector{Vector{Int64}} # Indices of elemental K in global S.nzval
    n::Int64 #number of DOF total
    ndofe::Int64 #number of DOF in elemental stiffness matrix
    cp::Vector{Int64} #S.colptr
    rv::Vector{Int64} #S.rowval
    nnz::Int64 #length(S.nzval)

    function TrussOptParams(nodeids, 
            dofids, 
            P, 
            freeids, 
            inzs, 
            n, 
            ndofe, 
            cp, 
            rv, 
            nnz)

        return new(nodeids, dofids, P, freeids, inzs, n, ndofe, cp, rv, nnz)
    end
end

function TrussOptParams(model::TrussModel)
    nodeids = getproperty.(model.elements, :nodeIDs)
    dofids = getproperty.(model.elements, :globalID)
    P = model.P
    freeids = model.freeDOFs
    inzs = allinz(model)
    n = model.nDOFs
    ndofe = 6
    cp = model.S.colptr
    rv = model.S.rowval
    nnz = length(model.S.nzval)

    return TrussOptParams(nodeids, dofids, P, freeids, inzs, n, ndofe, cp, rv, nnz)
end

mutable struct OptTrussNode
    node::TrussNode
    x::Float64
    y::Float64
    z::Float64
    activity::Vector{Bool}
end

function OptTrussNode(node::TrussNode, activity::Vector{Bool})
    @assert length(activity) == 3

    return OptTrussNode(node, node.position..., activity)
end

nodeActivityDict = Dict(
    :x => 1,
    :X => 1,
    :y => 2,
    :Y => 2, 
    :z => 3,
    :Z => 3)

function OptTrussNode(node::TrussNode, activity::Vector{Symbol})
    @assert length(activity) <= 3
    @assert all([in(sym, keys(nodeActivityDict)) for sym in activity])

    active = [false, false, false]
    for sym in activity
        active[nodeActivityDict[sym]] = true
    end

    return OptTrussNode(node, node.position..., active)

end

mutable struct OptTrussElement
    element::TrussElement
    A::Float64
    E::Float64
    activity::Vector{Bool}
end

function OptTrussElement(element::TrussElement, activity::Vector{Bool})
    @assert length(activity) == 2

    return OptTrussElement(element, element.section.A, element.section.E, activity)
end