export MfuncDict
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

export VfuncDict
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

export DfuncDict
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

export release2DOF
release2DOF = Dict(:fixedfixed => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    :freefixed => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    :fixedfree => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    :freefree => [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    :joist => [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0])

export planarDOFs
planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

