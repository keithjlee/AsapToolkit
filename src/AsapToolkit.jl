"""
    AsapToolkit

Companion utilities for [Asap.jl](https://github.com/keithjlee/Asap):

- `SteelSections` — the AISC v15 shapes database (`W`, `C`, `L`, `LL`, `WT`,
  `HSSRect`, `HSSRound`) with converters to Asap sections
  (`toASAPframe`, `toASAPtruss`)
- `AsapSections` — polygonal cross-section geometry and analysis
  (areas, centroids, depth maps)
- IO — JSON topology export (`topologize`) and Grasshopper interop (`GHsave`)

The parametric structure generators (`Warren2D`, `SpaceFrame`, `Frame`, …),
plot-geometry extraction (`Geo`, `ElementDisplacements`), FDM translation
(`to_network`), and model utilities (`clear_supports!`,
`element_connectivity`) moved into Asap itself for v1.0 — load Asap to use
them.
"""
module AsapToolkit
using Reexport
using Asap, LinearAlgebra
using JSON

include("SteelSections/SteelSections.jl")

include("IO/topologize.jl")
include("IO/GH/GHAsap.jl")
export GHsave

include("AsapSections/AsapSections.jl")
@reexport using .AsapSections

end # module AsapToolkit
