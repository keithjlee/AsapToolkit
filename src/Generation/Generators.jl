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

Optional inputs:
- `diaphragm::Bool = false` Include floor diaphragm (as X braces between columns)
- `diaphragmSection::Section = nothing` Element section for diaphragm members (defaults to braceSection)
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
        baseRelease = :pinned,
        diaphragm = false,
        diaphragmSection = nothing,
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
        if last(node.position) > last(base)
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
    #diaphragm
    if diaphragm

        dsection = isnothing(diaphragmSection) ? braceSection : diaphragmSection

        for k = 2:nz+1
            for i = 1:nx
                for j = 1:ny
                    diaph1 = Element(nodes[i, j, k], nodes[i+1, j+1, k], dsection, :freefree)
                    diaph2 = Element(nodes[i+1, j, k], nodes[i, j+1, k], dsection, :freefree)

                    diaph1.id = diaph2.id = :diaphragm
                    push!(braces, diaph1, diaph2)
                end
            end
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

function generatewarren2d(xpositions::Vector{<:Real},
    ypositions::Vector{<:Real},
    ypositions2::Vector{<:Real},
    section::Asap.AbstractSection;
    type = :arch,
    base = [0., 0., 0.])

    @assert length(xpositions) == length(ypositions) == length(ypositions2) + 1
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
    else
        longid = :bottomchord
        shortid = :topchord
    end

    n = length(xpositions)
    i = 1

    ## generate long chord up to symmetry
    for (x, y) in zip(xpositions, ypositions)

        node = TrussNode([x, y, 0.] .+ base, :free)
        if i == 1
            node.dof = [false, false, false]
            node.id = :pin
            i += 1
        else
            node.id = longid
        end

        push!(nodes, node)
        push!(longids, count)
        count += 1
    end

    ## generate other side of symmetry
    Lhalf = last(xpositions)
    base2 = [Lhalf, 0., 0.] .+ base
    incs = Lhalf .- reverse(xpositions[1:end-1])
    for (inc, y) in zip(incs, reverse(ypositions[1:end-1]))
        node = TrussNode(base2 .+ [inc, y, 0.], :free)
        node.id = longid
        push!(nodes, node)
        push!(longids, count)
        count += 1
    end
    fixnode!(last(nodes), :yfixed)
    last(nodes).id = :roller


    #generate shorter chord
    shortids = Vector{Int64}()
    for i = 1:length(ypositions2)

        xposition = mean(xpositions[i:i+1])

        node = TrussNode([xposition, ypositions2[i], 0.] .+ base, :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)
        count += 1
    end

    #other side
    for (i, y) in zip(length(xpositions):length(longids), reverse(ypositions2))
        x = first(mean(getproperty.(nodes, :position)[i:i+1]))

        node = TrussNode([x, y, 0.], :free)
        node.id = shortid

        push!(nodes, node)
        push!(shortids, count)

        count += 1
    end


    #elements
    elements = Vector{TrussElement}()

    #long chords
    for i = 1:length(longids) - 1
        element = TrussElement(nodes, longids[i:i+1], section)
        element.id = longid

        push!(elements, element)
    end

    #short chords
    for i = 1:length(shortids) - 1
        element = TrussElement(nodes, shortids[i:i+1], section)
        element.id = shortid

        push!(elements, element)
    end

    #webs
    for i = 1:length(longids) - 1
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
    # truss = Warren2D(model, n, dx, dy, section, type, base)

    #output
    return model
end

struct SpaceFrame
    truss::TrussModel
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    dz::Real
    section::Asap.AbstractSection
    support::Symbol
    load::Vector{<:Real}
    base::Vector{<:Real}
    ibottom::Matrix{Int64}
    itop::Matrix{Int64}
    isquares::Matrix{Vector{Int64}}
    isupport::Vector{Int64}
    iX1::Vector{Int64}
    iX2::Vector{Int64}
    iY1::Vector{Int64}
    iY2::Vector{Int64}
end


function generatespaceframe(nx::Integer,
        dx::Real,
        ny::Integer,
        dy::Real,
        dz::Real,
        section::Asap.AbstractSection;
        support = :corner,
        load = [0., 0., -10.],
        base = [0., 0., 0.])

    #generate nodes for bottom plane
    bottomnodes = [TrussNode([dx * (i-1), dy * (j-1), 0.] .+ base, :free) for i in 1:nx+1, j in 1:ny+1]
    for node in bottomnodes
        node.id = :bottom
    end

    #generate top nodes
    xinit = dx / 2
    yinit = dy / 2

    topnodes = [TrussNode([dx * (i-1) + xinit, dy * (j-1) + yinit, dz], :free) for i in 1:nx, j in 1:ny]
    for node in topnodes
        node.id = :top
    end

    #elements
    elements = Vector{TrussElement}()

    #generate bottom horizontal elements
    #parallel to x
    for j = 1:ny+1
        for i = 1:nx
            element = TrussElement(bottomnodes[i,j], bottomnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx+1
        for j = 1:ny
            element = TrussElement(bottomnodes[i,j], bottomnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate top horizontal elements
    #parallel to x
    for j = 1:ny
        for i = 1:nx-1
            element = TrussElement(topnodes[i,j], topnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx
        for j = 1:ny-1
            element = TrussElement(topnodes[i,j], topnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate web elements
    for i = 1:nx
        for j = 1:ny
            e1 = TrussElement(topnodes[i,j], bottomnodes[i,j], section)
            e2 = TrussElement(topnodes[i,j], bottomnodes[i+1,j], section)
            e3 = TrussElement(topnodes[i,j], bottomnodes[i,j+1], section)
            e4 = TrussElement(topnodes[i,j], bottomnodes[i+1,j+1], section)

            e1.id = e2.id = e3.id = e4.id = :web

            push!(elements, e1, e2, e3, e4)
        end
    end


    #generate node index matrices
    count = 1

    ibottomnodes = zeros(Int64, nx+1, ny+1)
    for j in 1:ny+1
        for i in 1:nx+1
            ibottomnodes[i,j] = count
            count += 1
        end
    end

    itopnodes = zeros(Int64, nx, ny)
    for j in 1:ny
        for i in 1:nx
            itopnodes[i,j] = count
            count += 1
        end
    end

    isquares = [[ibottomnodes[i,j],
        ibottomnodes[i+1,j],
        ibottomnodes[i,j+1],
        ibottomnodes[i+1,j+1]] for i in 1:nx, j in 1:ny]


    #generate supports
    if support == :corner
        for i in [1, nx+1]
            for j in [1, ny+1]
                fixnode!(bottomnodes[i,j], :pinned)
                bottomnodes[i,j].id = :support
            end
        end
    elseif support == :center
        isupport = nx % 2 == 1 ? [Int(ceil(nx / 2))] : Int.([ceil(nx / 2), floor(nx /2)])
        jsupport = ny % 2 == 1 ? [Int(ceil(ny / 2))] : Int.([ceil(ny / 2), floor(ny /2)])

        iset = vcat([vec(isquares[i,j]) for i in isupport, j in jsupport]...)

        for i in iset
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    end

    flatnodes = [vec(bottomnodes); vec(topnodes)]

    #generate loads
    loads  = [NodeForce(node, load) for node in flatnodes if node.id != :support]

    #assemble
    truss = TrussModel(flatnodes, elements, loads)
    solve!(truss)

    isupport = findall(truss.nodes, :support)

    ix1 = ibottomnodes[1,:]
    ix2 = ibottomnodes[end,:]
    iy1 = ibottomnodes[:,1]
    iy2 = ibottomnodes[:,end]


    spaceframe = SpaceFrame(truss,
        nx,
        dx,
        ny,
        dy,
        dz,
        section,
        support,
        load,
        base,
        ibottomnodes,
        itopnodes,
        isquares,
        isupport,
        ix1,
        ix2,
        iy1,
        iy2)

    return spaceframe
end

function generatespaceframe(nx::Integer,
    dx::Real,
    ny::Integer,
    dy::Real,
    z0::Real,
    interpolator::Interpolations.AbstractExtrapolation,
    section::Asap.AbstractSection;
    support = :corner,
    load = [0., 0., -10.],
    base = [0., 0., 0.])

    @assert bounds(iterpolator.itp) == ((0.0, 1.0), (0.0, 1.0)) "Interpolator must be parameterized from 0 → 1 for both x,y coordinates"

    #generate nodes for bottom plane
    bottomnodes = [TrussNode([dx * (i-1), dy * (j-1), 0.] .+ base, :free) for i in 1:nx+1, j in 1:ny+1]
    for node in bottomnodes
        node.id = :bottom
    end

    #generate top nodes
    xinit = dx / 2
    yinit = dy / 2

    xmax = dx * nx
    ymax = dy * ny

    topnodes = [TrussNode([dx * (i-1) + xinit, 
        dy * (j-1) + yinit, 
        z0 + interpolator(dx * (i-1) / xmax, dy * (j-1) /ymax)], 
        :free) for i in 1:nx, j in 1:ny]

    for node in topnodes
        node.id = :top
    end

    #elements
    elements = Vector{TrussElement}()

    #generate bottom horizontal elements
    #parallel to x
    for j = 1:ny+1
        for i = 1:nx
            element = TrussElement(bottomnodes[i,j], bottomnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx+1
        for j = 1:ny
            element = TrussElement(bottomnodes[i,j], bottomnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate top horizontal elements
    #parallel to x
    for j = 1:ny
        for i = 1:nx-1
            element = TrussElement(topnodes[i,j], topnodes[i+1,j], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #parallel to y
    for i = 1:nx
        for j = 1:ny-1
            element = TrussElement(topnodes[i,j], topnodes[i,j+1], section)
            element.id = :bottom

            push!(elements, element)
        end
    end

    #generate web elements
    for i = 1:nx
        for j = 1:ny
            e1 = TrussElement(topnodes[i,j], bottomnodes[i,j], section)
            e2 = TrussElement(topnodes[i,j], bottomnodes[i+1,j], section)
            e3 = TrussElement(topnodes[i,j], bottomnodes[i,j+1], section)
            e4 = TrussElement(topnodes[i,j], bottomnodes[i+1,j+1], section)

            e1.id = e2.id = e3.id = e4.id = :web

            push!(elements, e1, e2, e3, e4)
        end
    end


    #generate node index matrices
    count = 1

    ibottomnodes = zeros(Int64, nx+1, ny+1)
    for j in 1:ny+1
        for i in 1:nx+1
            ibottomnodes[i,j] = count
            count += 1
        end
    end

    itopnodes = zeros(Int64, nx, ny)
    for j in 1:ny
        for i in 1:nx
            itopnodes[i,j] = count
            count += 1
        end
    end

    isquares = [[ibottomnodes[i,j],
        ibottomnodes[i+1,j],
        ibottomnodes[i,j+1],
        ibottomnodes[i+1,j+1]] for i in 1:nx, j in 1:ny]


    #generate supports
    if support == :corner
        for i in [1, nx+1]
            for j in [1, ny+1]
                fixnode!(bottomnodes[i,j], :pinned)
                bottomnodes[i,j].id = :support
            end
        end
    elseif support == :center
        isupport = nx % 2 == 1 ? [Int(ceil(nx / 2))] : Int.([ceil(nx / 2), floor(nx /2)])
        jsupport = ny % 2 == 1 ? [Int(ceil(ny / 2))] : Int.([ceil(ny / 2), floor(ny /2)])

        iset = vcat([vec(isquares[i,j]) for i in isupport, j in jsupport]...)

        for i in iset
            bottomnodes[i].id = :support
            fixnode!(bottomnodes[i], :pinned)
        end
    end

    flatnodes = [vec(bottomnodes); vec(topnodes)]

    #generate loads
    loads  = [NodeForce(node, load) for node in flatnodes if node.id != :support]

    #assemble
    truss = TrussModel(flatnodes, elements, loads)
    solve!(truss)

    isupport = findall(truss.nodes, :support)


    ix1 = ibottomnodes[1,:]
    ix2 = ibottomnodes[end,:]
    iy1 = ibottomnodes[:,1]
    iy2 = ibottomnodes[:,end]


    spaceframe = SpaceFrame(truss,
        nx,
        dx,
        ny,
        dy,
        dz,
        section,
        support,
        load,
        base,
        ibottomnodes,
        itopnodes,
        isquares,
        isupport,
        ix1,
        ix2,
        iy1,
        iy2)

    return spaceframe
end