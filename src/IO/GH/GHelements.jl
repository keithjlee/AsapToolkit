struct GHsection
    E::Float64
    G::Float64
    A::Float64
    Ix::Float64
    Iy::Float64
    J::Float64

    function GHsection(section::TrussSection)

        return new(section.E, 1., section.A, 1., 1., 1.)

    end

    function GHsection(section::Section)

        return new(section.E, section.G, section.A, section.Ix, section.Iy, section.J)

    end
end

const release2bool = Dict(
    :fixedfixed => [false, false, false, false, false, false],
    :freefixed => [true, true, true, false, false, false],
    :fixedfree => [false, false, false, true, true, true],
    :freefree => [true, true, true, true, true, true],
    :joist => [false, true, true, false, true, true]
)

struct GHelement
    iStart::Int64
    iEnd::Int64
    elementID::Int64
    section::GHsection
    release::Vector{Bool}
    psi::Float64
    localx::Vector{Float64}
    localy::Vector{Float64}
    localz::Vector{Float64}
    forces::Vector{Float64}
    axialforce::Float64
    id::String

    function GHelement(element::TrussElement)
        istart, iend = element.nodeIDs .- 1
        elementID = element.elementID - 1
        section = GHsection(element.section)
        psi = element.Ψ
        lx, ly, lz = element.LCS
        id = isnothing(element.id) ? "" : string(element.id)
        forces = element.forces
        axialforce = forces[2]
        release = fill(true, 6)

        return new(istart, iend, elementID, section, release, psi, lx, ly, lz, forces, axialforce, id)
    end
    
    function GHelement(element::Element)
        istart, iend = element.nodeIDs .- 1
        elementID = element.elementID - 1
        section = GHsection(element.section)
        psi = element.Ψ
        lx, ly, lz = element.LCS
        id = isnothing(element.id) ? "" : string(element.id)
        forces = element.forces
        axialforce = forces[7]
        release = release2bool[element.release]

        return new(istart, iend, elementID, section, release, psi, lx, ly, lz, forces, axialforce, id)
    end
end
