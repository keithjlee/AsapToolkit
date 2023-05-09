function Imin(P::T, L::T, E::T; ϕ = 0.85) where T <: Real
    ϕ * P * L^2 / pi^2 / E
end

function Amin(P::T, σ::T; ϕ = 0.85) where T <: Real
    ϕ * P / σ
end

function trusssizer(element::TrussElement, σ::Real)
    F = element.forces[2]
    L = element.length
    E = element.section.E

    if F < 0
        return [-Amin(F, σ), -Imin(F, L, E)]
    else
        return [Amin(F, σ), 0.]
    end

end