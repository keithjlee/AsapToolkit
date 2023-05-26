function localvector(x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)
    [x2 - x1, y2 - y1, z2 - z1]
end

function localvector(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, id::Vector{Int64})
    [X[id[2]] - X[id[1]], Y[id[2]] - Y[id[1]], Z[id[2]] - Z[id[1]]]
end


"""
Local stiffness matrix for truss element
"""
function klocal(E::Float64, A::Float64, L::Float64)
    E * A / L * [1 -1; -1 1]
end

"""
Length of element
"""
function L(x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end

"""
Transformation matrix of 3d truss element
"""
function Rtruss(Cx::Float64, Cy::Float64, Cz::Float64)
    [Cx Cy Cz 0. 0. 0.; 0. 0. 0. Cx Cy Cz]
end

"""
elemental stiffness matrix in GCS
"""
function kglobal(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Float64, A::Float64, id::Vector{Int64})

    i1, i2 = id

    veclocal = [X[i2] - X[i1], Y[i2] - Y[i1], Z[i2] - Z[i1]]
    len = norm(veclocal)

    cx, cy, cz = veclocal ./ len
    r = Rtruss(cx, cy, cz)
    kloc = klocal(E, A, len)

    r' * kloc * r
end

"""
Assemble the global stiffness matrix from a vector of elemental stiffness matrices
"""
function assembleglobalK(elementalKs::Vector{Matrix{Float64}}, p::TrussOptParams)

    nz = zeros(p.nnz)

    for (k, i) in zip(elementalKs, p.inzs)
        nz[i] .+= vec(k)
    end

    SparseMatrixCSC(p.n, p.n, p.cp, p.rv, nz)
end

"""
Solve for the displacements of the FREE DOFs of the system
"""
function solveU(K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end

"""
    replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

Replace the values of `values[indices]` with the values in `newvalues` in a differentiable way. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function replacevalues(values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})
    
    v2 = copy(values)
    v2[indices] .= newvalues

    return v2
end

"""
    addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

Add the values of `increments` to the current values in `values` at `indices`. Does NOT perform any bounds checking or vector length consistency. This should be done before calling this function.
"""
function addvalues(values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    v2 = copy(values)
    v2[indices] .+= increments

    return v2
end