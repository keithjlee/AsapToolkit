abstract type AbstractParams end

mutable struct OptParams <: AbstractParams
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

    function OptParams(nodeids, 
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

function OptParams(model::TrussModel)
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

    return OptParams(nodeids, dofids, P, freeids, inzs, n, ndofe, cp, rv, nnz)
end

mutable struct OptTrussNode
    varid::Int64
    x::Float64
    y::Float64
    z::Float64
    activity::Vector{Bool}
end

mutable struct OptTrussElement
    varid::Int64
    E::Float64
    A::Float64
    L::Float64
    activity::Vector{Bool}
end
