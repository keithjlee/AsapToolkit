using XLSX
include("utilities.jl")

include("types.jl")
export W, C, L, LL, WT, HSSRect, HSSRound

include("functions.jl")
export allW, allC, allL, allLL, allWT, allHSSRect, allHSSRound

export Wnames, Cnames, Lnames, LLnames, WTnames, HSSRectnames, HSSRoundnames

export toASAPframe, toASAPtruss