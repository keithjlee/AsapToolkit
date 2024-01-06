module AsapSections

using LinearAlgebra

include("Sections.jl")
export SolidSection, VoidSection, CompoundSection, OffsetSection

include("SectionProperties.jl")
export SectionProperties

include("GeneralFunctions.jl")
export poly_area
export center_at_centroid!
export rotate_section!
export translate_section!
export move!

include("DepthAnalysis.jl")
export sutherland_hodgman
export sutherland_hodgman_abs
export intersection
export depth_map
export area_from_depth
export depth_from_area

end # module AsapSections
