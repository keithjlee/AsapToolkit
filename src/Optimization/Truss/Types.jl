abstract type AbstractOptParams end
abstract type AbstractVariable end
abstract type TrussOptVariable end
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

mutable struct SpatialVariable <: TrussOptVariable
    i::Int64 #index of node, e.g. X[i] is the spatial variable
    val::Float64 #value
    lb::Float64 #lower bound of variable
    ub::Float64 #upper bound of variable
    axis::Symbol #which spatial coordinate?
    iglobal::Int64 # position in the vector of active design variables

    function SpatialVariable(nodeindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        return new(nodeindex, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::TrussNode, value::Float64, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end

    function SpatialVariable(node::TrussNode, lowerbound::Float64, upperbound::Float64, axis::Symbol = :Z)

        @assert in(axis, validaxes)

        value = node.position[axis2ind[axis]]

        return new(node.nodeID, value, lowerbound, upperbound, axis)
    end
end

mutable struct AreaVariable <: TrussOptVariable
    i::Int64 #index of element, e.g. A[i] is the area variable
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Int64

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
    referencevariable::TrussOptVariable

    function CoupledVariable(node::TrussNode, ref::SpatialVariable)
        return new(node.nodeID, ref)
    end

    function CoupledVariable(element::TrussElement, ref::AreaVariable)
        return new(element.elementID, ref)
    end
end

const TrussVariables = Union{SpatialVariable, AreaVariable, CoupledVariable}

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


function populate!(indexer::TrussOptIndexer, var::SpatialVariable)
    field_local, field_global = axis2field[var.axis]

    push!(getfield(indexer, field_local), var.i)
    push!(getfield(indexer, field_global), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::AreaVariable)
    push!(getfield(indexer, :iA), var.i)
    push!(getfield(indexer, :iAg), var.iglobal)

end

function populate!(indexer::TrussOptIndexer, var::CoupledVariable)
    if typeof(var.referencevariable) == SpatialVariable
        field_local, field_global = axis2field[var.referencevariable.axis]

        push!(getfield(indexer, field_local), var.i)
        push!(getfield(indexer, field_global), var.referencevariable.iglobal)
    else
        push!(getfield(indexer, :iA), var.i)
        push!(getfield(indexer, :iAg), var.referencevariable.iglobal)
    end
end

# create the active variable → truss variable indexer
function TrussOptIndexer(vars::Vector{TrussVariables})
    indexer = TrussOptIndexer(Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Int64}())

    for var in vars
        populate!(indexer, var)
    end
    
    return indexer
end

mutable struct TrussOptProblem <: AbstractOptProblem
    model::TrussModel #the reference truss model for optimization
    params::TrussOptParams #optimization params
    indexer::TrussOptIndexer #pointers to design variables and full variables
    variables::Vector{TrussVariables}
    X::Vector{Float64} #all X coordinates |n_node|
    Y::Vector{Float64} #all Y coordinates |n_node|
    Z::Vector{Float64} #all Z coordinates |n_node|
    E::Vector{Float64} #all element young's modulii |n_element|
    A::Vector{Float64} #all element areas |n_element|
    values::Vector{Float64} #design variables
    lb::Vector{Float64} #lower bounds of variables
    ub::Vector{Float64} #upper bounds of variables

    function TrussOptProblem(model::TrussModel, variables::Vector{TrussVariables})
        @assert model.processed

        #extract global parameters
        xyz = Asap.nodePositions(model)
        X = xyz[:, 1]; Y = xyz[:, 2]; Z = xyz[:, 3]
        E = getproperty.(getproperty.(model.elements, :section), :E)
        A = getproperty.(getproperty.(model.elements, :section), :A)

        #extract secondary parameters
        params = TrussOptParams(model)

        #assign global id to variables
        vals = Vector{Float64}()
        lowerbounds = Vector{Float64}()
        upperbounds = Vector{Float64}()

        #assign an index to all unique variables, collect value and bounds
        i = 1
        for var in variables
            if typeof(var) <: TrussOptVariable
                var.iglobal  = i
                i += 1
                push!(vals, var.val)
                push!(lowerbounds, var.lb)
                push!(upperbounds, var.ub)
            end
        end

        #generate indexer between design variables and truss parameters
        indexer = TrussOptIndexer(variables)

        #generate a truss optimization problem
        return new(model, 
            params, 
            indexer, 
            variables, 
            X, 
            Y, 
            Z, 
            E, 
            A, 
            vals, 
            lowerbounds, 
            upperbounds)

    end
end
