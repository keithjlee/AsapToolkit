"""
displacement function for the transverse translation and in-plane rotation for a GIVEN PLANE IN THE LCS OF AN ELEMENT:

u_xy= [ v₁
        θ₁
        v₂
        θ₂ ]

IE for the local XY plane:

v₁, v₂ are the start and end displacements in the local Y direction
θ₁, θ₂ are the start and end rotations in the **local Z** direction (ie rotation in plane of local XY)

Gives:

v_y(x) = N × u_xy (translational displacement in local Y at point x)

"""
function N(x::Float64, L::Float64)
    n1 = 1 - 3(x/L)^2 + 2(x/L)^3
    n2 = x * (1 - x/L)^2
    n3 = 3(x/L)^2 - 2(x/L)^3
    n4 = x^2/L * (-1 + x/L)

    return [n1 n2 n3 n4]
end

"""
Axial displacement function: linear interpolation between start and end displacements
"""
function Naxial(x::Float64, L::Float64)
    n1 = 1 - x/L
    n2 = x / L

    return [n1 n2]
end

"""
    displacements(element::Element; n::Integer = 20)

Get the [3 × n] matrix where each column represents the local [x,y,z] displacement of the element from end forces
"""
function unodal(element::Element; n::Integer = 20)

    # base properties
    ulocal = element.R * [element.nodeStart.displacement; element.nodeEnd.displacement] .* release2DOF[element.release]
    L = element.length

    # extracting relevant nodal DOFs
    uX = ulocal[[1, 7]]
    uY = ulocal[[2, 6, 8, 12]]
    uZ = ulocal[[3, 5, 9, 11]] .* [1, -1, 1, -1]

    # discretizing length of element
    xrange = range(0, L, n)

    nA = vcat(Naxial.(xrange, L)...)
    nT = vcat(N.(xrange, L)...)

    dx = nA * uX
    dy = nT * uY
    dz = nT * uZ

    # [dx' ; dy' ; dz']
    dx, dy, dz
end

"""
Accumlate the internal forces cause by a line load to an element
"""
function accumulatedisp!(
    load::LineLoad, 
    xvals::Vector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    E = load.element.section.E
    Istrong = load.element.section.Izz
    Iweak = load.element.section.Iyy

    release = load.element.release

    # distributed load magnitudes in LCS
    wx, wy, wz = (R  * load.value) .* [1, -1, -1]

    # extract relevant function
    dfunction = DLineLoad[release]

    Dy .-= dfunction.(wy, L, xvals, E, Istrong)
    Dz .-= dfunction.(wz, L, xvals, E, Iweak)
end

"""
Accumlate the internal forces cause by a point load to an element
"""
function accumulatedisp!(
    load::PointLoad, 
    xvals::Vector{Float64}, 
    Dy::Vector{Float64},
    Dz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    E = load.element.section.E
    Istrong = load.element.section.Izz
    Iweak = load.element.section.Iyy
    release = load.element.release
    frac = load.position

    # distributed load magnitudes in LCS
    px, py, pz = (R  * load.value) .* [1, -1, -1]

    # extract relevant function
    dfunction = DPointLoad[release]

    Dy .-= dfunction.(py, L, xvals, frac, E, Istrong)
    Dz .-= dfunction.(pz, L, xvals, frac, E, Iweak)
end

"""
    ulocal(element::Element, model::Model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in LCS
"""
function ulocal(element::Element, model::Model; resolution = 20)

    L = element.length

    xinc = collect(range(0, L, resolution))

    D = unodal(element; n = resolution)

    for load in model.loads[element.loadIDs]
        accumulatedisp!(load, xinc, D[2,:], D[3,:])
    end

    return D
end

"""
    uglobal(element::Element, model::Model; resolution = 20)

Get the [3 × resolution] matrix of xyz displacements in GCS
"""
function uglobal(element::Element, model::Model; resolution = 20)

    L = element.length

    xinc = collect(range(0, L, resolution))

    D = unodal(element; n = resolution)

    for load in model.loads[element.loadIDs]
        accumulatedisp!(load, xinc, D[2,:], D[3,:])
    end

    return hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)
end

struct ElementDisplacements
    element::Element
    resolution::Integer
    x::Vector{Float64}
    ulocal::Matrix{Float64}
    uglobal::Matrix{Float64}
    basepositions::Matrix{Float64}
end

"""
    ElementDisplacements(element::Element, model::Model; resolution = 20)

Get the local/global displacements of an element
"""
function ElementDisplacements(element::Element, model::Model; resolution = 20)
    L = element.length

    xinc = collect(range(0, L, resolution))

    Dx, Dy, Dz = unodal(element; n = resolution)

    for load in model.loads[element.loadIDs]
        accumulatedisp!(load, xinc, Dy, Dz)
    end

    D = [Dx'; Dy'; Dz']

    Dglobal = hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)

    basepoints = element.nodeStart.position .+ first(element.LCS) * xinc'

    return ElementDisplacements(element, resolution, xinc, D, Dglobal, basepoints)
end

function ElementDisplacements(elements::Vector{<:Asap.FrameElement}, model::Model; resolution = 20)

    xstore = Vector{Float64}()
    ulocalstore = Vector{Matrix{Float64}}()
    uglobalstore = Vector{Matrix{Float64}}()
    basepointstore = Vector{Matrix{Float64}}()

    resolution = Int(round(resolution / length(elements)))

    for element in elements
        L = element.length

        xinc = collect(range(0, L, resolution))

        Dx, Dy, Dz = unodal(element; n = resolution)

        for load in model.loads[element.loadIDs]
            accumulatedisp!(load, xinc, Dy, Dz)
        end

        D = [Dx'; Dy'; Dz']
        Dglobal = hcat([sum(Δ .* element.LCS) for Δ in eachcol(D)]...)

        basepoints = element.nodeStart.position .+ first(element.LCS) * xinc'

        if isempty(xstore)
            xstore = [xstore; xinc]
        else
            xstore = [xstore; xstore[end] .+ xinc]
        end

        push!(ulocalstore, D)
        push!(uglobalstore, Dglobal)
        push!(basepointstore, basepoints)
    end

    return ElementDisplacements(elements[1], resolution, xstore, hcat(ulocalstore...), hcat(uglobalstore...), hcat(basepointstore...))
end

"""
    displacements(model::Model, increment::Real)

Get the displacements of all elements in a model
"""
function displacements(model::Model, increment::Real)
    results = Vector{ElementDisplacements}()

    ids = groupbyid(model.elements)

    for id in ids
        elements = model.elements[id]
        L = sum(getproperty.(elements, :length))
        n = max(Int(round(L/increment)), 2)

        push!(results, ElementDisplacements(elements, model))
    end

    return results
end