abstract type AbstractOptParams end
abstract type AbstractVariable end
abstract type AbstractOptProblem end
abstract type AbstractIndexer end

struct TrussOptParams <: AbstractOptParams
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

mutable struct SpatialVariable <: AbstractVariable
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

mutable struct AreaVariable <: AbstractVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64

    function AreaVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(elementindex, value, lowerbound, upperbound)
    end

    function AreaVariable(element::TrussElement, value::Float64, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, value, lowerbound, upperbound)
    end

    function AreaVariable(element::TrussElement, lowerbound::Float64, upperbound::Float64)
        return new(element.elementID, element.section.A, lowerbound, upperbound)
    end
end

mutable struct CoupledVariable <: AbstractVariable
    i::Int64
    referencevariable::AbstractVariable

    function CoupledVariable(node::TrussNode, ref::SpatialVariable)
        return new(node.nodeID, ref)
    end

    function CoupledVariable(element::TrussElement, ref::AreaVariable)
        return new(element.elementID, ref)
    end
end

mutable struct TrussOptIndexer <: AbstractIndexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
end

# quick reference to relevant field in TrussOptIndexer from variables
axis2field = Dict(:X => (:iX, :iXg),
    :x => (:iX, :iXg),
    :Y => (:iY, :iYg),
    :y => (:iY, :iYg),
    :Z => (:iZ, :iZg),
    :z => (:iZ, :iZg))

function populate!(indexer::TrussOptIndexer, var::SpatialVariable, i::Int64)
    field_local, field_global = axis2field(var.axis)

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), i)
end

function populate!(indexer::TrussOptIndexer, var::AreaVariable, i::Int64)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), i)
end

function populate!(indexer::TrussOptIndexer, var::CoupledVariable, i::Int64)
    
end

function TrussOptIndexer(vars::Vector{Union{SpatialVariable, AreaVariable}})
    indexer = TrussOptIndexer(Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}())

    for (i, var) in enumerate(vars)
        populate!(indexer, var, i)
    end
    
    return indexer
end

mutable struct TrussOptProblem <: AbstractOptProblem
    model::TrussModel #the reference truss model for optimization
    params::TrussOptParams #optimization params
    variables::Vector{Union{SpatialVariable, AreaVariable}}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    values::Vector{Float64} #design variables
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables
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
