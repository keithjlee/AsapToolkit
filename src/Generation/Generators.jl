abstract type AbstractGenerator end

include("Frame.jl"); export Frame
include("Spaceframe.jl"); export SpaceFrame
include("Warren.jl"); export Warren2D
include("SpaceframeBeam.jl"); export SpaceFrameBeam
include("BakerTruss.jl"); export BakerTruss
include("TrussFrame.jl"); export TrussFrame
include("GridNetwork.jl"); export GridNetwork