export MfuncDict
const MfuncDict = Dict((LineLoad, :fixedfixed) => MLine_fixedfixed,
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

const MLineLoad = Dict(:fixedfixed => MLine_fixedfixed,
    :freefree => MLine_freefree,
    :fixedfree => MLine_fixedfree,
    :freefixed => MLine_freefixed,
    :joist => MLine_freefree
    )

const MPointLoad = Dict(:fixedfixed => MPoint_fixedfixed,
    :freefree => MPoint_freefree,
    :fixedfree => MPoint_fixedfree,
    :freefixed => MPoint_freefixed,
    :joist => MPoint_freefree)

export VfuncDict
const VfuncDict = Dict((LineLoad, :fixedfixed) => VLine_fixedfixed,
    (LineLoad, :freefree) => VLine_freefree,
    (LineLoad, :fixedfree) => VLine_fixedfree,
    (LineLoad, :freefixed) => VLine_freefixed,
    (LineLoad, :joist) => VLine_freefree,
    (PointLoad, :fixedfixed) => VPoint_fixedfixed,
    (PointLoad, :freefree) => VPoint_freefree,
    (PointLoad, :fixedfree) => VPoint_fixedfree,
    (PointLoad, :freefixed) => VPoint_freefixed,
    (PointLoad, :joist) => VPoint_freefree)

const VLineLoad = Dict(:fixedfixed => VLine_fixedfixed,
    :freefree => VLine_freefree,
    :fixedfree => VLine_fixedfree,
    :freefixed => VLine_freefixed,
    :joist => VLine_freefree)

const VPointLoad = Dict(:fixedfixed => VPoint_fixedfixed,
    :freefree => VPoint_freefree,
    :fixedfree => VPoint_fixedfree,
    :freefixed => VPoint_freefixed,
    :joist => VPoint_freefree)

export DfuncDict
const DfuncDict = Dict((LineLoad, :fixedfixed) => DLine_fixedfixed,
    (LineLoad, :freefree) => DLine_freefree,
    (LineLoad, :fixedfree) => DLine_fixedfree,
    (LineLoad, :freefixed) => DLine_freefixed,
    (LineLoad, :joist) => DLine_freefree,
    (PointLoad, :fixedfixed) => DPoint_fixedfixed,
    (PointLoad, :freefree) => DPoint_freefree,
    (PointLoad, :fixedfree) => DPoint_fixedfree,
    (PointLoad, :freefixed) => DPoint_freefixed,
    (PointLoad, :joist) => DPoint_freefree)

const DLineLoad = Dict(:fixedfixed => DLine_fixedfixed,
    :freefree => DLine_freefree,
    :fixedfree => DLine_fixedfree,
    :freefixed => DLine_freefixed,
    :joist => DLine_freefree)

const DPointLoad = Dict(:fixedfixed => DPoint_fixedfixed,
    :freefree => DPoint_freefree,
    :fixedfree => DPoint_fixedfree,
    :freefixed => DPoint_freefixed,
    :joist => DPoint_freefree)
    
export release2DOF
const release2DOF = Dict(:fixedfixed => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    :freefixed => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    :fixedfree => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    :freefree => [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    :joist => [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0])

export planarDOFs
const planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

"""
Collect all elements that are part of the same continuous member (ie collect shattered elements)
"""
function groupbyid(elements::Vector{<:Asap.AbstractElement})
    ichecked = Vector{Int64}()
    indices = Vector{Vector{Int64}}()

    elementids = getproperty.(elements, :elementID)

    for id in elementids
        
        in(id, ichecked) && continue

        igroup = findall(elementids .== id)

        push!(indices, igroup)
        push!(ichecked, id)
    end

    return indices
end