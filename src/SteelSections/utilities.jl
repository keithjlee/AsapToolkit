data = XLSX.readdata(joinpath(@__DIR__, "data/aisc-shapes-database-v15.0.xlsx"), "Database v15.0!CH2:ED2041")

Wrange = 1:283
Crange = 352:383
Lrange = 424:560
WTrange = 561:843
LLrange = 886:1524
HSSRectrange = 1525:1912
HSSRoundrange = 1913:2040
names = data[:,1]

# column index dictionary
colDict = Dict(
    :name => 1,
    :W => 2,
    :A => 3,
    :d => 4,
    :b => 12,
    :bf => 9,
    :tw => 14,
    :tf => 17,
    :t => 19,
    :Ix => 36,
    :Zx => 37,
    :Sx => 38,
    :rx => 39,
    :Iy => 40,
    :Zy => 41,
    :Sy => 42,
    :ry => 43,
    :J => 47,
    :Cw => 48,
    :Ht => 6,
    :B => 11,
    :tHSS => 21,
    :OD => 8
    )

Wfields = [:name,
    :A,
    :d,
    :bf,
    :tw,
    :tf,
    :Ix,
    :Zx,
    :Sx,
    :rx,
    :Iy,
    :Zy,
    :Sy,
    :ry,
    :J,
    :Cw]

Wfactors = [1., 
    1., 
    1., 
    1., 
    1., 
    1e6, 
    1e3, 
    1e3, 
    1., 
    1e6, 
    1e3, 
    1e3, 
    1., 
    1e3, 
    1e9]


Cfactors = [1., 1., 1., 1., 1., 1e6, 1e3, 1e3, 1., 1e6, 1e3, 1e3, 1., 1e3, 1e9]


Lfields = [:name,
    :A,
    :b,
    :d,
    :t,
    :Ix,
    :Zx,
    :Sx,
    :rx,
    :Iy,
    :Zy,
    :Sy,
    :ry,
    :J,
    :Cw]

Lfactors = [1., 1., 1., 1., 1e6, 1e3, 1e3, 1., 1e6, 1e3, 1e3, 1., 1e3, 1e9]

LLfields = [:name,
    :A,
    :b,
    :d,
    :t,
    :Ix,
    :Zx,
    :Sx,
    :rx,
    :Iy,
    :Zy,
    :Sy,
    :ry]

LLfactors = [1, 1, 1, 1, 1e6, 1e3, 1e3, 1, 1e6, 1e3, 1e3, 1]

HSSRectfields = [:name,
    :A,
    :Ht,
    :B,
    :tHSS,
    :Ix,
    :Zx,
    :Sx,
    :rx,
    :Iy,
    :Zy,
    :Sy,
    :ry,
    :J]

HSSRectfactors = [1, 1, 1, 1, 1e6, 1e3, 1e3, 1, 1e6, 1e3, 1e3, 1, 1e3]

HSSRoundfields = [:name,
    :A,
    :OD,
    :tHSS,
    :Ix,
    :Zx,
    :Sx,
    :rx,
    :Iy,
    :Zy,
    :Sy,
    :ry,
    :J]

HSSRoundfactors = [1, 1, 1, 1e6, 1e3, 1e3, 1, 1e6, 1e3, 1e3, 1, 1e3]