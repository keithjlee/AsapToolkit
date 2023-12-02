abstract type AbstractGeo end

include("TrussGeometry.jl"); export TrussGeo
include("ModelGeometry.jl"); export ModelGeo
include("NetworkGeometry.jl"); export NetworkGeo

export Geo
function Geo end

Geo(model::TrussModel) = TrussGeo(model)
Geo(model::Model) = ModelGeo(model)
Geo(network::Network) = NetworkGeo(network)