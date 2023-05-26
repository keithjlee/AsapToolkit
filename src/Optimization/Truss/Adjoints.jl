function ChainRulesCore.rrule(::typeof(localvector), x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)

    vector = localvector(x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)

    function localvector_pullback(v̄)
        return (NoTangent(),
            dot(v̄, [-1., 0., 0.]),
            dot(v̄, [1., 0., 0.]),
            dot(v̄, [0., -1., 0.]),
            dot(v̄, [0., 1., 0.]),
            dot(v̄, [0., 0., -1.]),
            dot(v̄, [0., 0., 1.]))
    end

    return vector, localvector_pullback
end

"""
Adjoint w/r/t element variables E, A, L for the local stiffness matrix of a truss element
"""
function ChainRulesCore.rrule(::typeof(klocal), E::Float64, A::Float64, L::Float64)
    k = klocal(E, A, L)

    function klocal_pullback(k̄)

        ∇E = dot(k̄, (A / L * [1 -1; -1 1]))
        ∇A = dot(k̄, (E / L * [1 -1; -1 1]))
        ∇L = dot(k̄, (- E * A / L^2 * [1 -1; -1 1]))

        return (NoTangent(), ∇E, ∇A, ∇L)
    end

    return k, klocal_pullback
end

"""
Explicit length adjoint -- not necessary but included
"""
function ChainRulesCore.rrule(::typeof(L), x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)
    len = L(x1::Float64, x2::Float64, y1::Float64, y2::Float64, z1::Float64, z2::Float64)

    function L_pullback(dL)
        xdiff = x2 - x1
        ydiff = y2 - y1
        zdiff = z2 - z1

        return (NoTangent(), -dL * xdiff / len, 
            dL * xdiff / len,
            -dL * ydiff / len,
            dL * ydiff / len,
            -dL * zdiff / len,
            dL * zdiff / len)
    end

    return len, L_pullback
end

"""
Adjoint for global transformation matrix
"""
function ChainRulesCore.rrule(::typeof(Rtruss), Cx::Float64, Cy::Float64, Cz::Float64)
    R = Rtruss(Cx, Cy, Cz)

    function Rtruss_pullback(R̄)
        return (NoTangent(),
            dot(R̄, [1. 0. 0. 0. 0. 0.; 0. 0. 0. 1. 0. 0.]),
            dot(R̄, [0. 1. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 0.]),
            dot(R̄, [0. 0. 1. 0. 0. 0.; 0. 0. 0. 0. 0. 1.]))
    end

    return R, Rtruss_pullback
end

"""
The sensitivity of K w/r/t an elemental Kₑ is the proportional stiffness added to K from Kₑ

Output is a vector: [nElements × [nDOFe × nDOFe]] of elemental stiffness matrix sensitivites
"""
function ChainRulesCore.rrule(::typeof(assembleglobalK), Eks::Vector{Matrix{Float64}}, p::TrussOptParams)
    K = assembleglobalK(Eks, p)

    function K_pullback(K̄)
        dEks = [Matrix(K̄[id, id]) for id in p.dofids]

        return NoTangent(), dEks, NoTangent()
    end

    return K, K_pullback
end

"""
u = inv(K) * P

if obj = f(u), then the gradient of obj with respect to an independent variable x is achieved through the chain rule:

dObj/dx = df/du ⋅ du/dK ⋅ ... = ū ⋅ dy/dK ⋅ ...

For this rule, we are concerned with finding du/dK, or the [ndof × ndof] matrix of sensitivites that we can propagate backwards to the final objective.

Given df/du = [ndof × 1] = ū is the gradient of the objective function with respect to displacements u, the sensitivity is:

du/dK = - uᵀ ⊗ K⁻¹
df/dK = du/dK ū = - (uᵀ ⊗ K⁻¹)ū

Can be rearranged such that:
ΔK = K⁻¹ū

df/dK = -uᵀ ⊗ ΔK

Which is an [ndof × ndof] matrix where:

Columnᵢ = uᵢ .* ΔK
"""
function ChainRulesCore.rrule(::typeof(solveU), K::SparseMatrixCSC{Float64, Int64}, p::TrussOptParams)
    u = solveU(K, p)

    function solveU_pullback(ū)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], ū)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solveU_pullback
end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement.

df/dnewvalues = df/dvalues ⋅ dvalues / dnewvalues = v̄ ⋅ dvalues / dnewvalues

Where v̄ = [nvalues × 1], dvalues/dnewvalues = [nvalues × nnewvalues] so:

df/dnewvalues = (dvalues/dnewvalues)ᵀv̄ = [nnewvalues × 1]

Is simply the values of v̄ at the indices of the new values.
"""
function ChainRulesCore.rrule(::typeof(replacevalues), values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

    v = replacevalues(values, indices, newvalues)

    function replacevalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, replacevalues_pullback 
end


"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement
"""
function ChainRulesCore.rrule(::typeof(addvalues), values::Vector{Float64}, indices::Vector{Int64}, increments::Vector{Float64})

    v = addvalues(values, indices, increments)

    function addvalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, addvalues_pullback 
end