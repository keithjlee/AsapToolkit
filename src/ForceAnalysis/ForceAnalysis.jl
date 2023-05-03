#### LINE - free free
export MLine_freefree
function MLine_freefree(w::T, L::T, x::T) where T <: Real
    w * x / 2 * (L - x)
end

export VLine_freefree
function VLine_freefree(w::T, L::T, x::T) where T <: Real
    w * (L / 2 - x)
end

export DLine_freefree
function DLine_freefree(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w * x / 24 / E / I * (L^3 - 2L*x^2 + x^3)
end

#### LINE - fixed fixed
export MLine_fixedfixed
function MLine_fixedfixed(w::T, L::T, x::T) where T <: Real
    w / 12 * (6 * L * x - L^2 - 6 * x^2)
end

export VLine_fixedfixed
function VLine_fixedfixed(w::T, L::T, x::T) where T <: Real
    w * (L / 2 - x)
end

export DLine_fixedfixed
function DLine_fixedfixed(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w * x ^ 2 / 24 / E / I * (L - x)^2
end


#### LINE - free fixed
export MLine_freefixed
function MLine_freefixed(w::T, L::T, x::T) where T <: Real
    R1 = 3 * w * L / 8

    R1 * x - w * x^2 / 2
end

export VLine_freefixed
function VLine_freefixed(w::T, L::T, x::T) where T <: Real
    R1 = 3 * w * L / 8

    R1 - w * x
end

export DLine_freefixed
function DLine_freefixed(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w * x / 48 / E / I * (L^3 - 3L * x^2 + 2x^3)
end

### LINE - fixed free

export MLine_fixedfree
function MLine_fixedfree(w::T, L::T, x::T) where T <: Real
    - MLine_freefixed(w, L, L-x)
end

export VLine_fixedfree
function VLine_fixedfree(w::T, L::T, x::T) where T <: Real
    - VLine_freefixed(w, L, L-x)
end

export DLine_fixedfree
function DLine_fixedfree(w::T, L::T, x::T, E::T, I::T) where T <: Real
    - DLine_freefixed(w, L, L-x, E, I)
end

#### POINT - freefree
export MPoint_freefree
function MPoint_freefree(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    x < a ? P * b * x / L : P * b * a / L - P * a * (x - a) / L
end

export VPoint_freefree
function VPoint_freefree(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    x < a ? P * b / L : -P * a / L
end

export DPoint_freefree
function DPoint_freefree(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P * b * x / 6 / E / I / L * (L^2 - b^2 - x^2)
    else
        P * a * (L - x) / 6 / E / I / L * (2L * x - x^2 - a^2)
    end
end

export MPoint_fixedfixed
function MPoint_fixedfixed(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    R1 = P * b^2 / L^3 * (3a + b)
    R2 = P * a^2 / L^3 * (a + 3b)

    if x < a
        R1 * x - P * a * b^2 / L^2
    else
        2 * P * a^2 * b^2 / L^3 - R2 * (x - a)
    end
end

export VPoint_fixedfixed
function VPoint_fixedfixed(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    x < a ? P * b^2 / L^3 * (3a + b) : - P * a^2 / L^3 * (a + 3b)
end

export DPoint_fixedfixed
function DPoint_fixedfixed(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P * b^2 * x^2 / 6 / E / I / L^3 * (3a * L - 3a * x - b * x)
    else
        P * a^3 * b^3 / 3 / E / I / L^3 - P * a^2 * (b + a - x)^2 / (6E * I * L^3) * (3b*L -3b*(b+a-x) - a*(b+a-x))
    end

end