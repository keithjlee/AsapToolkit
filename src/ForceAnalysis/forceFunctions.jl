function PLine(w::T, L::T, x::T) where T <: Real
    @assert 0 ≤ x ≤ L

    R = - w * L / 2
    # assumes start reaction is in line with local x
    -(R + w * x)
end

"""
    MLine_freefree(w, L, x)

Simply supported internal moment at `x` given distributed load `w` on beam of length `L`
"""
function MLine_freefree(w, L, x)
    @assert 0 ≤ x ≤ L

    w * x / 2 * (L - x)
end

MLine(::Element{Asap.FreeFree}, w, L, x) = MLine_freefree(w, L, x)
MLine(::Element{Asap.Joist}, w, L, x) = MLine_freefree(w, L, x)

"""
    VLine_freefree(w, L, x)

Simply supported internal shear at `x` given distributed load `w` on beam of length `L`
"""
function VLine_freefree(w, L, x)
    @assert 0 ≤ x ≤ L

    w * (L / 2 - x)
end

VLine(::Element{Asap.FreeFree}, w, L, x) = VLine_freefree(w, L, x)
VLine(::Element{Asap.Joist}, w, L, x) = VLine_freefree(w, L, x)

"""
    DLine_freefree(w, L, x, E, I)

Simply supported transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DLine_freefree(w, L, x, E, I)
    @assert 0 ≤ x ≤ L

    w * x / 24 / E / I * (L^3 - 2L*x^2 + x^3)
end

DLine(::Element{Asap.FreeFree}, w, L, x, E, I) = DLine_freefree(w, L, x, E, I)
DLine(::Element{Asap.Joist}, w, L, x, E, I) = DLine_freefree(w, L, x, E, I)


"""
    ThetaLine_freefree(w, L, x, E, I)

Simply supported rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaLine_freefree(w, L, x, E, I)
    @assert 0 ≤ x ≤ L

    w/24/E/I * (L^3 - 6L * x^2 + 4x^3)
end

ThetaLine(::Element{Asap.FreeFree}, w, L, x, E, I) = ThetaLine_freefree(w, L, x, E, I)
ThetaLine(::Element{Asap.Joist}, w, L, x, E, I) = ThetaLine_freefree(w, L, x, E, I)

"""
    MLine_fixedfixed(w, L::T, x::T) where T <: Real

Fixed-fixed internal moment at `x` given distributed load `w` on beam of length `L`
"""
function MLine_fixedfixed(w, L, x)
    @assert 0 ≤ x ≤ L

    w / 12 * (6 * L * x - L^2 - 6 * x^2)
end

MLine(::Element{Asap.FixedFixed}, w, L, x) = MLine_fixedfixed(w, L, x)


"""
    VLine_fixedfixed(w, L, x)

Simply supported internal shear at `x` given distributed load `w` on beam of length `L`
"""
function VLine_fixedfixed(w, L, x)
    @assert 0 ≤ x ≤ L

    w * (L / 2 - x)
end

VLine(::Element{Asap.FixedFixed}, w, L, x) = VLine_fixedfixed(w, L, x)

"""
    DLine_fixedfixed(w, L, x, E, I)

Simply supported transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DLine_fixedfixed(w, L, x, E, I)
    @assert 0 ≤ x ≤ L

    w * x ^ 2 / 24 / E / I * (L - x)^2
end

DLine(::Element{Asap.FixedFixed}, w, L, x, E, I) = DLine_fixedfixed(w, L, x, E, I)

"""
    ThetaLine_fixedfixed(w, L, x, E, I)

Simply supported rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaLine_fixedfixed(w, L, x, E, I)
    @assert 0 ≤ x ≤ L

    w/24/E/I * (L^2 * x^2 - 2L*x^3 +x^4)
end

ThetaLine(::Element{Asap.FixedFixed}, w, L, x, E, I) = ThetaLine_fixedfixed(w, L, x, E, I)

"""
    MLine_freefixed(w, L, x)

Free-fixed internal moment at `x` given distributed load `w` on beam of length `L`
"""
function MLine_freefixed(w, L, x)
    @assert 0 ≤ x ≤ L

    R1 = 3 * w * L / 8

    R1 * x - w * x^2 / 2
end

MLine(::Element{Asap.FreeFixed}, w, L, x) = MLine_freefixed(w, L, x)


"""
    VLine_freefixed(w, L, x)

Free-fixed internal shear at `x` given distributed load `w` on beam of length `L`
"""
function VLine_freefixed(w, L, x)
    @assert 0 ≤ x ≤ L

    R1 = 3 * w * L / 8

    R1 - w * x
end

VLine(::Element{Asap.FreeFixed}, w, L, x) = VLine_freefixed(w, L, x)

"""
    DLine_freefixed(w, L, x, E, I)

Free-fixed transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DLine_freefixed(w, L, x, E, I)
    @assert 0 ≤ x ≤ L
    w * x / 48 / E / I * (L^3 - 3L * x^2 + 2x^3)
end

DLine(::Element{Asap.FreeFixed}, w, L, x, E, I) = DLine_freefixed(w, L, x, E, I)


"""
    ThetaLine_freefixed(w, L, x, E, I)

Free-fixed rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaLine_freefixed(w, L, x, E, I)
    @assert 0 ≤ x ≤ L
    w/48/E/I * (L^3 - 9L*x^2 + 8x^3)
end

ThetaLine(::Element{Asap.FreeFixed}, w, L, x, E, I) = ThetaLine_freefixed(w, L, x, E, I)

"""
    MLine_fixedfree(w, L, x)

Fixed-free internal moment at `x` given distributed load `w` on beam of length `L`
"""
function MLine_fixedfree(w, L, x)
    MLine_freefixed(w, L, L-x)
end

MLine(::Element{Asap.FixedFree}, w, L, x) = MLine_fixedfree(w, L, x)

"""
    VLine_fixedfree(w, L, x)

Fixed-free internal shear at `x` given distributed load `w` on beam of length `L`
"""
function VLine_fixedfree(w, L, x) 
    -VLine_freefixed(w, L, L-x)
end

VLine(::Element{Asap.FixedFree}, w, L, x) = VLine_fixedfree(w, L, x)

"""
    DLine_fixedfree(w, L, x, E, I)

Fixed-free transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DLine_fixedfree(w, L, x, E, I)
    DLine_freefixed(w, L, L-x, E, I)
end

DLine(::Element{Asap.FixedFree}, w, L, x, E, I) = DLine_fixedfree(w, L, x, E, I)

"""
    ThetaLine_fixedfree(w, L, x, E, I)

Fixed-free rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaLine_fixedfree(w, L, x, E, I)
    -ThetaLine_freefixed(w, L, (L-x), E, I)
end

ThetaLine(::Element{Asap.FixedFree}, w, L, x, E, I) = ThetaLine_fixedfree(w, L, x, E, I)

#### POINT - freefree

function PPoint(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    R = - P / 2

    x < a ? -R : -(R + P)
end

PPoint(::Element, P, L, x, frac) = PPoint(P, L, x, frac)

"""
    MPoint_freefree(P, L, x, frac)

Simply supported internal moment at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function MPoint_freefree(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    b = L - a

    x < a ? P * b * x / L : P * b * a / L - P * a * (x - a) / L
end

MPoint(::Element{Asap.FreeFree}, P, L, x, frac) = MPoint_freefree(P, L, x, frac)

"""
    VPoint_freefree(P, L, x, frac)

Simply supported internal shear at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function VPoint_freefree(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    b = L - a

    x < a ? P * b / L : -P * a / L
end

VPoint(::Element{Asap.FreeFree}, P, L, x, frac) = VPoint_freefree(P, L, x, frac)

"""
    DPoint_freefree(P, L, x, frac, E, I)

Simply supported transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DPoint_freefree(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    b = L - a

    if x < a
        P * b * x / 6 / E / I / L * (L^2 - b^2 - x^2)
    else
        P * a * (L - x) / 6 / E / I / L * (2L * x - x^2 - a^2)
    end
end

DPoint(::Element{Asap.FreeFree}, P, L, x, frac, E, I) = DPoint_freefree(P, L, x, frac, E, I)

"""
    ThetaPoint_freefree(P, L, x, frac, E, I)

Simply supported rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaPoint_freefree(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L

    a = frac * L
    b = L - a

    if x < a
        P * b / 6 / E / I / L * (L^2 - b^2 - 3x^2)
    else
        P * a / 6 / E / I / L * (3x^2 - 6L*x +2L^2 +a^2)
    end
end

ThetaPoint(::Element{Asap.FreeFree}, P, L, x, frac, E, I) = ThetaPoint_freefree(P, L, x, frac, E, I)

"""
    MPoint_fixedfixed(P, L, x, frac)

Fixed-fixed internal moment at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function MPoint_fixedfixed(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
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

MPoint(::Element{Asap.FixedFixed}, P, L, x, frac) = MPoint_fixedfixed(P, L, x, frac)


"""
    VPoint_fixedfixed(P, L, x, frac)

Fixed-fixed internal shear at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function VPoint_fixedfixed(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    b = L - a

    x < a ? P * b^2 / L^3 * (3a + b) : - P * a^2 / L^3 * (a + 3b)
end

VPoint(::Element{Asap.FixedFixed}, P, L, x, frac) = VPoint_fixedfixed(P, L, x, frac)


"""
    DPoint_fixedfixed(P, L, x, frac, E, I)


Fixed-fixed transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DPoint_fixedfixed(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    
    a = frac * L
    b = L - a

    if x < a
        P * b^2 * x^2 / 6 / E / I / L^3 * (3a * L - 3a * x - b * x)
    else
        P * a^2 * (L - x)^2 / (6E * I * L^3) * (3b*L -3b*(L-x) - a*(L-x))
    end

end

DPoint(::Element{Asap.FixedFixed}, P, L, x, frac, E, I) = DPoint_fixedfixed(P, L, x, frac, E, I)

"""
    ThetaPoint_fixedfixed(P, L, x, frac, E, I)


Fixed-fixed rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaPoint_fixedfixed(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L

    a = frac * L
    b = L - a

    if x < a
        P * b^2 / 6 / E / I / L^3 * (6a * L * x - x^2 * (9a + 3b))
    else
        P * a^2 / (6E * I * L^3) * (3b*L^2 +a*L^2 - 6b*L^2 + (L-2x) * (6b*L + 2a*L) +6b*L*x - (2L*x - 3x^2) * (3b + a))
    end

end

ThetaPoint(::Element{Asap.FixedFixed}, P, L, x, frac, E, I) = ThetaPoint_fixedfixed(P, L, x, frac, E, I)

"""
    MPoint_freefixed(P, L, x, frac)

Free-fixed internal moment at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function MPoint_freefixed(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L

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

MPoint(::Element{Asap.FreeFixed}, P, L, x, frac) = MPoint_freefixed(P, L, x, frac)

"""
    MPoint_fixedfree(P, L, x, frac)

Fixed-free internal moment at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function MPoint_fixedfree(P, L, x, frac)
    MPoint_freefixed(P, L, (L-x), 1-frac)
end

MPoint(::Element{Asap.FixedFree}, P, L, x, frac) = MPoint_fixedfree(P, L, x, frac)

"""
    VPoint_freefixed(P, L, x, frac)

Free-fixed internal shear at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function VPoint_freefixed(P, L, x, frac)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L

    a = frac * L
    b = L - a

    R1 = P*b^2/2/L^3 * (a + 2L)
    R2 = P*a/2/L^3 * (3L^2 - a^2)

    x < a ? R1 : -R2
end

VPoint(::Element{Asap.FreeFixed}, P, L, x, frac) = VPoint_freefixed(P, L, x, frac)

"""
    VPoint_fixedfree(P, L, x, frac)

Fixed-free internal shear at `x` given point load `P` on beam length `L` at a point `frac × L`
"""
function VPoint_fixedfree(P, L, x, frac)
    - VPoint_freefixed(P, L, (L-x), 1-frac)
end

VPoint(::Element{Asap.FixedFree}, P, L, x, frac) = VPoint_fixedfree(P, L, x, frac)

"""
    DPoint_freefixed(P, L, x, frac, E, I)

Free-fixed transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DPoint_freefixed(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    @assert 0 ≤ x ≤ L
    a = frac * L
    b = L - a

    if x < a
        P * b^2 * x / (12*E*I*L^3) * (3a * L^2 - 2L*x^2 -a*x^2)
    else
        P * a / (12*E*I*L^3) * (L-x)^2 * (3L^2*x - a^2*x - 2a^2*L)
    end

end

DPoint(::Element{Asap.FreeFixed}, P, L, x, frac, E, I) = DPoint_freefixed(P, L, x, frac, E, I)

"""
    DPoint_fixedfree(P, L, x, frac, E, I)

Fixed-free transverse displacement at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function DPoint_fixedfree(P, L, x, frac, E, I)
    DPoint_freefixed(P, L, (L-x), 1-frac, E, I)
end

DPoint(::Element{Asap.FixedFree}, P, L, x, frac, E, I) = DPoint_fixedfree(P, L, x, frac, E, I)

"""
    ThetaPoint_freefixed(P, L, x, frac, E, I)

Free-fixed rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaPoint_freefixed(P, L, x, frac, E, I)
    @assert 0 ≤ frac ≤ 1.
    a = frac * L
    b = L - a

    if x < a
        P*b^2/(12*E*I*L^3) * (3a*L^2 - x^2 * (6L + 3a))
    else
        P*a/(12*E*I*L^3) * (x^2*(9L^2-3a^2) - 12L^3*x + 4a^2*L^2 + 3L^4 - L^2*a^2)
    end

end

ThetaPoint(::Element{Asap.FreeFixed}, P, L, x, frac, E, I) = ThetaPoint_freefixed(P, L, x, frac, E, I)

"""
    ThetaPoint_fixedfree(P, L, x, frac, E, I)


Fixed-free rotation at `x` given:
- Distributed load `w`
- Beam length `L`
- Young's Modulus `E`
- Moment of Inertia `I`
"""
function ThetaPoint_fixedfree(P, L, x, frac, E, I)
    -ThetaPoint_freefixed(P, L, (L-x), 1-frac, E, I)
end

ThetaPoint(::Element{Asap.FixedFree}, P, L, x, frac, E, I) = ThetaPoint_fixedfree(P, L, x, frac, E, I)