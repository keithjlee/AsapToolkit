#### LINE - free free
export MLine_freefree
function MLine_freefree(w, L::T, x::T) where T <: Real
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

export ThetaLine_freefree
function ThetaLine_freefree(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w/24/E/I * (L^3 - 6L * x^2 + 4x^3)
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

export ThetaLine_fixedfixed
function ThetaLine_fixedfixed(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w/24/E/I * (L^2 * x^2 - 2L*x^3 +x^4)
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

export ThetaLine_freefixed
function ThetaLine_freefixed(w::T, L::T, x::T, E::T, I::T) where T <: Real
    w/48/E/I * (L^3 - 9L*x^2 + 8x^3)
end

### LINE - fixed free

export MLine_fixedfree
function MLine_fixedfree(w::T, L::T, x::T) where T <: Real
    -MLine_freefixed(w, L, L-x)
end

export VLine_fixedfree
function VLine_fixedfree(w::T, L::T, x::T) where T <: Real
    -VLine_freefixed(w, L, L-x)
end

export DLine_fixedfree
function DLine_fixedfree(w::T, L::T, x::T, E::T, I::T) where T <: Real
    DLine_freefixed(w, L, L-x, E, I)
end

export ThetaLine_fixedfree
function ThetaLine_fixedfree(w::T, L::T, x::T, E::T, I::T) where T <: Real
    # w/48/E/I * (L^3 * x * (L - x) - 3L * (L^2 -4L*x + 3x^2 - 6 * (L-x)^2))
    -ThetaLine_freefixed(w, L, (L-x), E, I)
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

export ThetaPoint_freefree
function ThetaPoint_freefree(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P * b / 6 / E / I / L * (L^2 - b^2 - 3x^2)
    else
        P * a / 6 / E / I / L * (3x^2 - 6L*x +2L^2 +a^2)
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
        P * a^2 * (L - x)^2 / (6E * I * L^3) * (3b*L -3b*(L-x) - a*(L-x))
    end

end

export ThetaPoint_fixedfixed
function ThetaPoint_fixedfixed(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P * b^2 / 6 / E / I / L^3 * (6a * L * x - x^2 * (9a + 3b))
    else
        P * a^2 / (6E * I * L^3) * (3b*L^2 +a*L^2 - 6b*L^2 + (L-2x) * (6b*L + 2a*L) +6b*L*x - (2L*x - 3x^2) * (3b + a))
    end

end

export MPoint_freefixed
function MPoint_freefixed(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    R1 = P*b^2/2/L^3 * (a + 2L)
    R2 = P*a/2/L^3 * (3L^2 - a^2)

    if x < a
        R1 * x
    else
        R1 * x - P * (x - a)
    end
end

export MPoint_fixedfree
function MPoint_fixedfree(P::T, L::T, x::T, frac::Float64) where T <: Real
    - MPoint_freefixed(P, L, (L-x), 1-frac)
end

export VPoint_freefixed
function VPoint_freefixed(P::T, L::T, x::T, frac::Float64) where T <: Real
    a = frac * L
    b = L - a

    R1 = P*b^2/2/L^3 * (a + 2L)
    R2 = P*a/2/L^3 * (3L^2 - a^2)

    x < a ? R1 : -R2
end

export VPoint_fixedfree
function VPoint_fixedfree(P::T, L::T, x::T, frac::Float64) where T <: Real
    - VPoint_freefixed(P, L, (L-x), 1-frac)
end

export DPoint_freefixed
function DPoint_freefixed(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P * b^2 * x / (12*E*I*L^3) * (3a * L^2 - 2L*x^2 -a*x^2)
    else
        P * a / (12*E*I*L^3) * (L-x)^2 * (3L^2*x - a^2*x - 2a^2*L)
    end

end

export DPoint_fixedfree
function DPoint_fixedfree(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    DPoint_freefixed(P, L, (L-x), 1-frac, E, I)
end

export ThetaPoint_freefixed
function ThetaPoint_freefixed(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    a = frac * L
    b = L - a

    if x < a
        P*b^2/(12*E*I*L^3) * (3a*L^2 - x^2 * (6L + 3a))
    else
        P*a/(12*E*I*L^3) * (x^2*(9L^2-3a^2) - 12L^3*x + 4a^2*L^2 + 3L^4 - L^2*a^2)
    end

end

export ThetaPoint_fixedfree
function ThetaPoint_fixedfree(P::T, L::T, x::T, frac::Float64, E::T, I::T) where T <: Real
    -ThetaPoint_freefixed(P, L, (L-x), 1-frac, E, I)
end