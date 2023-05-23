abstract type AbstractParams end
abstract type OptVariable end

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

mutable struct TrussOptNode
    node::TrussNode
    x::Float64
    y::Float64
    z::Float64
    activity::Vector{Bool}
end

function TrussOptNode(node::TrussNode, activity::Vector{Bool})
    @assert length(activity) == 3

    return TrussOptNode(node, node.position..., activity)
end

axis2ind = Dict(
    :x => 1,
    :X => 1,
    :y => 2,
    :Y => 2, 
    :z => 3,
    :Z => 3)

function TrussOptNode(node::TrussNode, activity::Vector{Symbol})
    @assert length(activity) <= 3
    @assert all([in(sym, keys(axis2ind)) for sym in activity])

    active = [false, false, false]
    for sym in activity
        active[axis2ind[sym]] = true
    end

    return TrussOptNode(node, node.position..., active)

end

mutable struct TrussOptElement
    element::TrussElement
    A::Float64
    E::Float64
    activity::Vector{Bool}
end

function TrussOptElement(element::TrussElement, activity::Vector{Bool})
    @assert length(activity) == 2

    return TrussOptElement(element, element.section.A, element.section.E, activity)
end

mutable struct TrussOptProblem
    model::TrussModel #the reference truss model for optimization
    variables::Vector{<:OptVariable}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    params::TrussOptParams #optimization params
    nspatial::Int64 #number of spatial variables (3 × n_node)
    ninternal::Int64 #number of internal variables (n_element)

    function TrussOptProblem(model::TrussModel, variables::Vector{<:OptVariable})
        @assert model.processed

        xyz = Asap.nodePositions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]

        E = [e.section.E for e in model.elements]
        A = [e.section.A for e in model.elements]

        params = TrussOptParams(model)

        nspatial = model.nNodes * 3
        ninternal = model.nElements * 3

        return new(model, variables, X, Y, Z, E, A, params, nspatial, ninternal)

    end
end


mutable struct SpatialVariable <: OptVariable
    i::Int64 #index of node, e.g. X[i] is the spatial variable
    val::Float64 #value
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    axis::Symbol #which spatial coordinate?

    function SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, keys(axis2ind))

        return new(nodeindex, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::TrussNode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, keys(axis2ind))

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::TrussNode, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, keys(axis2ind))

        value = node.position[axis2ind[axis]]

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end
end

mutable struct InternalVariable <: OptVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64

    function InternalVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(elementindex, value, lowerbound, upperbound)
    end

    function InternalVariable(element::TrussElement, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, value, lowerbound, upperbound)
    end

    function InternalVariable(element::TrussElement, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, element.section.A, lowerbound, upperbound)
    end
end

#generate objective function
mutable struct TrussOptIndexer
    iX::Vector{Int64}
    iY::Vector{Int64}
    iZ::Vector{Int64}
    iA::Vector{Int64}
    vX::Vector{Float64}
    vY::Vector{Float64}
    vZ::Vector{Float64}
    nXstart::Int64
    nYstart::Int64
    nZstart::Int64
    nEstart::Int64
    iVars::Vector{Int64}
    iVals::Vector{Float64}
    ubs::Vector{Float64}
    lbs::Vector{Float64}
end
