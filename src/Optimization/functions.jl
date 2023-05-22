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
function kglobal(posStart::Vector{Float64}, posEnd::Vector{Float64}, E::Float64, A::Float64)
    
    veclocal = posEnd .- posStart
    len = norm(veclocal)

    cx, cy, cz = veclocal ./ len
    r = Rtruss(cx, cy, cz)
    kloc = klocal(E, A, len)

    r' * kloc * r
end

function kglobal(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Float64, A::Float64, id::Vector{Int64})

    i1, i2 = id

    veclocal = [X[i2] - X[i1], Y[i2] - Y[i1], Z[i2] - Z[i1]]
    len = norm(veclocal)

    cx, cy, cz = veclocal ./ len
    r = Rtruss(cx, cy, cz)
    kloc = klocal(E, A, len)

    r' * kloc * r
end

function kglobal(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64}, E::Vector, A::Vector{Float64}, id::Vector{Int64})

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
function assembleglobalK(elementalKs::Vector{Matrix{Float64}}, p::AbstractParams)

    nz = zeros(p.nnz)

    for (k, i) in zip(elementalKs, p.inzs)
        nz[i] .+= vec(k)
    end

    SparseMatrixCSC(p.n, p.n, p.cp, p.rv, nz)
end


"""
Solve for the displacements of the FREE DOFs of the system
"""
function solveU(K::SparseMatrixCSC{Float64, Int64}, p::AbstractParams)
    id = p.freeids
    cg(K[id, id], p.P[id])
end