MfuncDict = Dict((LineLoad, :fixedfixed) => MLine_fixedfixed,
    (LineLoad, :freefree) => MLine_freefree,
    (LineLoad, :fixedfree) => MLine_fixedfree,
    (LineLoad, :freefixed) => MLine_freefixed,
    (LineLoad, :joist) => MLine_freefree,
    (PointLoad, :fixedfixed) => MPoint_fixedfixed,
    (PointLoad, :freefree) => MPoint_freefree,
    (PointLoad, :fixedfree) => MPoint_fixedfree,
    (PointLoad, :freefixed) => MPoint_freefixed,
    (PointLoad, :joist) => MPoint_freefree
    )

VfuncDict = Dict((LineLoad, :fixedfixed) => VLine_fixedfixed,
    (LineLoad, :freefree) => VLine_freefree,
    (LineLoad, :fixedfree) => VLine_fixedfree,
    (LineLoad, :freefixed) => VLine_freefixed,
    (LineLoad, :joist) => VLine_freefree,
    (PointLoad, :fixedfixed) => VPoint_fixedfixed,
    (PointLoad, :freefree) => VPoint_freefree,
    (PointLoad, :fixedfree) => VPoint_fixedfree,
    (PointLoad, :freefixed) => VPoint_freefixed,
    (PointLoad, :joist) => VPoint_freefree)

DfuncDict = Dict((LineLoad, :fixedfixed) => DLine_fixedfixed,
    (LineLoad, :freefree) => DLine_freefree,
    (LineLoad, :fixedfree) => DLine_fixedfree,
    (LineLoad, :freefixed) => DLine_freefixed,
    (LineLoad, :joist) => DLine_freefree,
    (PointLoad, :fixedfixed) => DPoint_fixedfixed,
    (PointLoad, :freefree) => DPoint_freefree,
    (PointLoad, :fixedfree) => DPoint_fixedfree,
    (PointLoad, :freefixed) => DPoint_freefixed,
    (PointLoad, :joist) => DPoint_freefree)

release2DOF = Dict(:fixedfixed => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    :freefixed => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    :fixedfree => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    :freefree => [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    :joist => [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0])

planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

"""
    N(x::Float64, L::Float64)

Cubic interpolation shape function value at position `x` on beam of length `L`. Multiply with end displacement vector [v₁, θ₁, v₂, θ₂]ᵀ, where v, θ are the transverse displacement and in-plane rotation for the XY or XZ plane
"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x * (1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    return [n1 n2 n3 n4]
end

"""
    N(x::Float64, L::Float64)

Linear interpolation shape function value at position `x` on beam of length `L`. Multiply with end displacement vector [u₁, u₂]ᵀ, where u is the axial displacement of the end nodes.
"""
function Naxial(x::Float64, L::Float64)
    n1 = 1 - x/L
    n2 = x / L

    return [n1 n2]
end

function dofdisplacement(element::Element; n::Integer = 10)
    
    # get the end node displacement vector in LCS
    u = element.R * ([element.nodeStart.displacement; element.nodeEnd.displacement] .* release2DOF[element.release])
    L = element.length

    xrange = range(0, L, n)

    uX = u[planarDOFs[:X]]
    uY = u[planarDOFs[:XY]]
    uZ = u[planarDOFs[:XZ]] #.* [1, -1, 1, -1]

    [@inbounds hcat([Naxial(x, L) * uX for x in xrange]...);
        @inbounds hcat([N(x, L) * uY for x in xrange]...);
        @inbounds hcat([N(x, L) * uZ for x in xrange]...)]


end