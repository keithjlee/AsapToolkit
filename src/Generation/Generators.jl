mutable struct Frame
    model::Model
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    nz::Integer
    dz::Real
    joistSpacing::Real
    columnSection::Section
    primarySection::Section
    joistSection::Section
    braceSection::Section
    columnRelease::Symbol
    primaryRelease::Symbol
    joistRelease::Symbol
    braceRelease::Symbol
    columnPsi::Real
    primaryPsi::Real
    joistPsi::Real
    base::Vector{<:Real}
    iExteriorXnodes::Vector{Integer}
    iExteriorYnodes::Vector{Integer}
    iExteriorXjoists::Vector{Integer}
    iExteriorYprimaries::Vector{Integer}
end

mutable struct Warren2D
    model::TrussModel
    n::Integer
    dx::Real
    dy::Real
    section::Asap.AbstractSection
    type::Symbol
    base::Vector{<:Real}
end

"""
    generateFrame(nx::Integer,...)

Generate a 3D frame model.

Required inputs:
- `nx::Integer` Number of bays in primary span direction
- `dx::Real` Primary span length
- `ny::Integer` Number of bays in secondary span direction
- `dy::Real` Secondary span length
- `nz::Integer` Number of stories
- `dz::Real` Story height
- `joistSpacing::Real` Approx. center spacing between secondary spans
- `columnSection::Section` Section of column elements
- `primarySection::Section` Section of primary beams
- `joistSection::Section` Section of secondary beams
- `braceSection::Section` Section for lateral bracing

Default inputs:
- `columnRelease::Symbol = :fixedfixed` DOF end release for column elements
- `primaryRelease::Symbol = :fixedfixed` DOF end release for primary elements
- `joistRelease::Symbol = :joist` DOF end release for joist elements
- `braceRelease::Symbol = :freefree` DOF end release for braces
- `columnPsi::Real = 0` Angle of roll Ψ for column LCS
- `primaryPsi::Real = π/2` Angle of roll Ψ for primary beam LCS
- `joistPsi::Real = π/2` Angle of roll Ψ for secondary beam LCS
- `base::Vector{Real} = [0, 0, 0]` base point for frame grid generation 
"""
function generateframe(nx::Integer,
        dx::Real,
        ny::Integer,
        dy::Real,
        nz::Integer,
        dz::Real,
        joistSpacing::Real,
        columnSection::Section,
        primarySection::Section,
        joistSection::Section,
        braceSection::Section;
        columnRelease = :fixedfixed,
        primaryRelease = :fixedfixed,
        joistRelease = :joist,
        braceRelease = :freefree,
        columnPsi = 0.,
        primaryPsi = pi/2,
        joistPsi = pi/2,
        base = [0., 0., 0.]
        )


    ########
    # nodes
    ########

    xoffset = [dx, 0., 0.]
    yoffset = [0., dy, 0.]
    zoffset = [0., 0., dz]

    #generate
    nodes = [Node(base + xoffset * i + yoffset * j + zoffset * k, :pinned) for i in 0:nx, j in 0:ny, k in 0:nz]

    #release non-ground nodes
    for node in nodes
        if last(node.position) > 0.
            fixnode!(node, :free)
        end
    end

    ######
    # columns
    #######
    #make columns
    columns = Vector{Element}()
    for i in 1:nx+1
        for j in 1:ny+1
            for k in 1:nz
                el = Element(nodes[i,j,k], nodes[i,j,k+1], columnSection, columnRelease)
                # el.Ψ = 0
                el.Ψ = columnPsi
                el.id = :column
                push!(columns, el)
            end
        end
    end


    #####
    # primary beams
    #####
    primaries = Vector{Element}()

    for k in 2:nz+1
        for i in 1:nx
            for j in 1:ny+1
                el = Element(nodes[i,j,k], nodes[i+1,j,k], primarySection, primaryRelease)
                el.id = :primary
                el.Ψ = primaryPsi
                push!(primaries, el)
            end
        end
    end

    ######
    # joists
    ######

    primreshaped = reshape(primaries, ny+1, nx, nz)

    secondaries = Vector{Union{BridgeElement, Element}}()
    njoists = Int(round(dx / joistSpacing))

    #main bridge elements
    for k in 1:nz
        for j = 1:nx
            for i = 1:ny
                bridges = [BridgeElement(primreshaped[i,j,k], x, primreshaped[i+1,j,k], x, joistSection, joistRelease) for x in range(0,1,njoists)[2:end-1]]

                for bridge in bridges
                    bridge.id = :joist
                    bridge.Ψ = joistPsi
                end

                push!(secondaries, bridges...)

            end
        end
    end

    #between column nodes
    for i in 1:nx+1
        for j in 1:ny
            for k in 2:nz+1
                bridge = Element(nodes[i,j,k], nodes[i,j+1,k], joistSection, joistRelease)
                # bridge.Ψ = 0.
                bridge.id = :joist
                bridge.Ψ = joistPsi
                push!(secondaries, bridge)
            end
        end
    end

    #########
    #braces
    #########

    braces = Vector{Element}()
    i = 1
    for j = [1, ny+1]
        for k = 1:nz
            brace1 = Element(nodes[i,j,k], nodes[i+1,j,k+1], braceSection, braceRelease)
            brace2 = Element(nodes[i+1,j,k], nodes[i,j,k+1], braceSection, braceRelease)

            brace1.id = brace2.id = :brace
            push!(braces, brace1, brace2)
        end
    end

    j = 1
    for i = [1, nx+1]
        for k = 1:nz
            brace1 = Element(nodes[i,j,k], nodes[i, j+1, k+1], braceSection, braceRelease)
            brace2 = Element(nodes[i,j+1,k], nodes[i,j,k+1], braceSection, braceRelease)

            brace1.id = brace2.id = :brace
            push!(braces, brace1, brace2)
        end
    end


    ######
    # dummy load and assembly
    ######

    loads = [LineLoad(j, [0., 0., -1.]) for j in secondaries]
    flatnodes = vec(nodes)
    elements = [columns; primaries; secondaries; braces]

    model = Model(flatnodes, elements, loads)
    solve!(model)


    #extract node/element indices
    iExteriorXnodes = [getproperty.(vec(nodes[:,1,:]), :nodeID); getproperty.(vec(nodes[:,end,:]), :nodeID)]
    iExteriorYnodes = [getproperty.(vec(nodes[1,:,:]), :nodeID); getproperty.(vec(nodes[end,:,:]), :nodeID)]

    iExteriorXjoists = Vector{Int64}()
    iExteriorYprimary = Vector{Int64}()

    xExtrema = [first(base), first(base) + dx * nx]
    yExtrema = [base[2], base[2] + dy * ny]
    for (i,element) in enumerate(model.elements)
        if element.id == :joist
            x = element.nodeStart.position[1]
            if minimum(abs.(xExtrema .- x)) <= model.tol
                push!(iExteriorXjoists, i)
            end
        elseif element.id == :primary
            y = element.nodeStart.position[2]
            if minimum(abs.(yExtrema .- y)) <= model.tol
                push!(iExteriorYprimary, i)
            end
        end
    end

    frame = Frame(model,
        nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        joistSpacing,
        columnSection,
        primarySection,
        joistSection,
        braceSection,
        columnRelease,
        primaryRelease,
        joistRelease,
        braceRelease,
        columnPsi,
        primaryPsi,
        joistPsi,
        base,
        iExteriorXnodes,
        iExteriorYnodes,
        iExteriorXjoists,
        iExteriorYprimary)

    return frame;
end

"""
    generatewarren2d(n::Integer,...)

Generate a 2D warren truss in the XY plane.

Required inputs:
- `n::Integer` Number of bays
- `dx::Real` Bay span
- `dy::Real` truss depth
- `section::Asap.AbstractSection` cross section of elements

Default inputs:
- `type::Symbol = :arch` :arch = long chord at bottom; :catenary = long chord at top
- `base::Vector{Real} = [0, 0, 0]` base point for truss generation
"""
function generatewarren2d(n::Integer,
        dx::Real,
        dy::Real,
        section::Asap.AbstractSection;
        type = :arch,
        base = [0., 0., 0.])

    @assert n % 2 != 0 "n must be odd"
    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{TrussNode}()

    #generate longer chord first
    if type == :arch
        longid = :topchord
        shortid = :bottomchord
        y = dy
    else
        longid = :bottomchord
        shortid = :topchord
        y = -dy
    end

    for i = 1:n
        xposition = dx * (i - 1)

        node = TrussNode([xposition, 0., 0.] .+ base, :free)
        if i == 1
            node.dof = [false, false, false]
            node.id = :pin
        elseif i == n
            node.dof = [true, false, false]
            node.id = :roller
        else
            node.id = longid
        end

        push!(nodes, node)
        push!(longids, count)

        count += 1
    end

    #generate shorter chord
    shortids = Vector{Int64}()
    x0 = dx / 2

    for i = 1:n-1
        xposition = x0 + dx * (i - 1)

        node = TrussNode([xposition, y, 0.] .+ base, :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)
        count += 1
    end

    #elements
    elements = Vector{TrussElement}()
    
    #long chords
    for i = 1:n-1
        element = TrussElement(nodes, longids[i:i+1], section)
        element.id = longid

        push!(elements, element)
    end

    #short chords
    for i = 1:n-2
        element = TrussElement(nodes, shortids[i:i+1], section)
        element.id = shortid

        push!(elements, element)
    end

    #webs
    for i = 1:n-1
        element = TrussElement(nodes, [longids[i], shortids[i]], section)
        element.id = :web
        push!(elements, element)

        element = TrussElement(nodes, [shortids[i], longids[i+1]], section)
        element.id = :web
        push!(elements, element)
    end

    #dummy load
    loads = [NodeForce(n, [0., -1., 0.],) for n in nodes[longid]]

    #assemble and solve
    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    #collect data
    truss = Warren2D(model, n, dx, dy, section, type, base)

    #output
    return truss
end