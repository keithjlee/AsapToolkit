struct SectionProperties
    section::PolygonalSection
    centroid::Vector{Float64}
    controid_c::Vector{Float64}
    centroid_t::Vector{Float64}
    area::Float64
    Ix::Float64
    Sx::Float64
    Zx::Float64
    rx::Float64
    Iy::Float64
    Sy::Float64
    Zy::Float64
    ry::Float64

    """
        SectionProperties(section::PolygonalSection)

    Get an extended section property analysis that includes:
    - centroid of compression/tension regions as defined by the horizontal plane passing through the global centroid
    - plastic modulii (Zx, Zy) about X and Y axes centered at the centroid
    - radii of gyration (rx, ry) about X and Y axes centered at the centroid
    """
    function SectionProperties(section::PolygonalSection)

        #radius of gyration
        rx = sqrt(section.Ix / section.area)
        ry = sqrt(section.Iy / section.area)

        #plastic modulus
        compression_section = sutherland_hodgman_abs(section, section.centroid[2]; return_section = true)

        Ac = compression_section.area
        At = section.area - Ac

        centroid_c = compression_section.centroid
        centroid_t = (section.area * section.centroid - compression_section.area * compression_section.centroid) / At

        ΔC = abs.(centroid_c .- section.centroid)
        ΔT = abs.(centroid_t .- section.centroid)

        Zx = Ac * ΔC[2] + At * ΔT[2]

        
        section2 = deepcopy(section)
        rotate_section!(section2, pi/2)

        #plastic modulus
        compression_section2 = sutherland_hodgman_abs(section2, section2.centroid[2]; return_section = true)

        Ac = compression_section2.area
        At = section2.area - Ac

        centroid_c = compression_section2.centroid
        centroid_t = (section2.area * section2.centroid - compression_section2.area * compression_section2.centroid) / At

        ΔC = abs.(centroid_c .- section2.centroid)
        ΔT = abs.(centroid_t .- section2.centroid)


        Zy = Ac * ΔC[2] + At * ΔT[2]

        return new(
            section,
            section.centroid,
            centroid_c,
            centroid_t,
            section.area,
            section.Ix,
            section.Sx,
            Zx,
            rx,
            section.Iy,
            section.Sy,
            Zy,
            ry
        )
    end

end
