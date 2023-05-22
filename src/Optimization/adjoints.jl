"""
Adjoint w/r/t element variables E, A, L for the local stiffness matrix of a truss element
"""
function ChainRulesCore.rrule(::typeof(klocal), E::Float64, A::Float64, L::Float64)
    k = klocal(E, A, L)

    function klocal_pullback(dk)

        ∇E = sum(dk .* (A / L * [1 -1; -1 1]))
        ∇A = sum(dk .* (E / L * [1 -1; -1 1]))
        ∇L = sum(dk .* (- E * A / L^2 * [1 -1; -1 1]))

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
        den = len^3

        xdiff = x2 - x1
        ydiff = y2 - y1
        zdiff = z2 - z1

        return (NoTangent(), -dL * xdiff / den, 
            dL * xdiff / den,
            -dL * ydiff / den,
            dL * ydiff / den,
            -dL * zdiff / den,
            dL * zdiff / den)
    end

    return len, L_pullback
end

"""
The sensitivity of K w/r/t an elemental Kₑ is the proportional stiffness added to K from Kₑ

Output is a vector: [nElements × [nDOFe × nDOFe]] of elemental stiffness matrix sensitivites
"""
function ChainRulesCore.rrule(::typeof(assembleglobalK), Eks::Vector{Matrix{Float64}}, p::AbstractParams)
    K = assembleglobalK(Eks, p)

    function K_pullback(dK)
        dEks = [Matrix(dK[id, id]) for id in p.dofids]

        return NoTangent(), dEks, NoTangent()
    end

    return K, K_pullback
end

"""
u = inv(K) * P

if obj = f(u), then the gradient of obj with respect to an independent variable x is achieved through the chain rule:

dObj/dx = df/du ⋅ du/dK ⋅ ...

For this rule, we are concerned with finding du/dK, or the [ndof × ndof] matrix of sensitivites that we can propagate backwards to the final objective.

Given df/du = [ndof × 1] = Δu is the gradient of the objective function with respect to displacements u, the sensitivity is:

du/dK = - uᵀ ⊗ K⁻¹
df/dK = du/dK ū = - (uᵀ ⊗ K⁻¹)Δu

Can be rearranged such that:
ΔK = K⁻¹Δu

df/dK = -uᵀ ⊗ ΔK

Which is an [ndof × ndof] matrix where:

Columnᵢ = uᵢ .* ΔK
"""
function ChainRulesCore.rrule(::typeof(solveU), K::SparseMatrixCSC{Float64, Int64}, p::AbstractParams)
    u = solveU(K, p)

    function solveU_pullback(Δ)

        #initialize
        dudK = zeros(p.n, p.n)

        #sensitivities w/r/t active DOFs
        dKΔ = cg(K[p.freeids, p.freeids], Δ)

        #assign to proper indices
        dudK[p.freeids, p.freeids] .= kron(u', dKΔ)

        return NoTangent(), -dudK, NoTangent()
    end

    return u, solveU_pullback
end

"""
Pullback of partial array replacement is simply the primal cotangent values *at* the indices of replacement
"""
function ChainRulesCore.rrule(::typeof(updatevalues), values::Vector{Float64}, indices::Vector{Int64}, newvalues::Vector{Float64})

    v = updatevalues(values, indices, newvalues)

    function updatevalues_pullback(v̄)

        return NoTangent(), NoTangent(), NoTangent(), v̄[indices]

    end

    return v, updatevalues_pullback 
end