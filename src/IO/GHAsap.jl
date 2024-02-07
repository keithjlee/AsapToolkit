struct GHnode
    position::Vector{Float64}
    dof::Vector{Bool}
    nodeID::Int64
    reaction::Vector{Float64}
    u::Vector{Float64}
    displacement::Vector{Float64}
    id::String

    function GHnode(node::TrussNode)
        position = node.position
        dof = [node.dof; [true, true, true]]
        nodeID = node.nodeID - 1
        reaction = node.reaction
        u = node.displacement
        displacement = node.displacement[1:3]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end

    function GHnode(node::Node)
        position = node.position
        dof = node.dof
        nodeID = node.nodeID - 1
        reaction = node.reaction
        u = node.displacement
        displacement = node.displacement[1:3]
        id = isnothing(node.id) ? "" : string(node.id)

        return new(position, dof, nodeID, reaction, u, displacement, id)
    end
end

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

struct GHelement
    iStart::Int64
    iEnd::Int64
    elementID::Int64
    section::GHsection
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
        id = isnothing(element.id) ? "" : element.id
        forces = element.forces
        axialforce = forces[2]

        return new(istart, iend, elementID, section, psi, lx, ly, lz, forces, axialforce, id)
    end
    
    function GHelement(element::Element)
        istart, iend = element.nodeIDs .- 1
        elementID = element.elementID - 1
        section = GHsection(element.section)
        psi = element.Ψ
        lx, ly, lz = element.LCS
        id = isnothing(element.id) ? "" : element.id
        forces = element.forces
        axialforce = forces[7]

        return new(istart, iend, elementID, section, psi, lx, ly, lz, forces, axialforce, id)
    end
end

abstract type GHload end
function GHload end

struct GHnodeforce <: GHload
    iNode::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::NodeForce)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : load.id

    return GHnodeforce(i, value, id)
end

struct GHnodemoment <: GHload
    iNode::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::NodeMoment)

    i = load.node.nodeID - 1
    value = load.value
    id = isnothing(load.id) ? "" : load.id

    return GHnodemoment(i, value, id)

end

struct GHlineload <: GHload
    iElement::Int64
    value::Vector{Float64}
    id::String
end

function GHload(load::LineLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : load.id

    return GHlineload(i, value, id)

end

struct GHpointload <: GHload
    iElement::Int64
    x::Float64
    value::Vector{Float64}
    id::String
end

function GHload(load::PointLoad)

    i = load.element.elementID - 1
    value = load.value
    id = isnothing(load.id) ? "" : load.id
    x = load.position

    return GHpointload(i, x, value, id)

end

struct GHmodel
    nodes::Vector{GHnode}
    elements::Vector{GHelement}
    loads::Vector{GHload}
    ifree::Vector{Int64}
    ifixed::Vector{Int64}