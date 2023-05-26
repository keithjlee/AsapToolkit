# valid indices for spatial variable assignment
const validaxes = [:X, :x, :Y, :y, :Z, :z]

# to indices in position vector
const axis2ind = Dict(:X => 1,
    :x => 1,
    :Y => 2,
    :y => 2,
    :Z => 3,
    :z => 3)

# quick reference to relevant field in TrussOptIndexer from variables
const axis2field = Dict(:X => (:iX, :iXg),
    :x => (:iX, :iXg),
    :Y => (:iY, :iYg),
    :y => (:iY, :iYg),
    :Z => (:iZ, :iZg),
    :z => (:iZ, :iZg))