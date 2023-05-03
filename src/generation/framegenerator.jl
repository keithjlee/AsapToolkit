function generateFrame(nx::Integer,
        dx::Real,
        ny::Integer,
        dy::Real,
        nz::Integer,
        dz::Real,
        joistspacing::Real,
        columnSection::Section,
        primarySection::Section,
        joistSection::Section,
        braceSection::Section;
        columnRelease = :fixedfixed,
        primaryRelease = :fixedfixed,
        joistRelease = :fixedfixed,
        braceRelease = :freefree,
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
                push!(primaries, el)
            end
        end
    end

    ######
    # joists
    ######

    primreshaped = reshape(primaries, ny+1, nx, nz)

    secondaries = Vector{Union{BridgeElement, Element}}()
    njoists = Int(round(dx / joistspacing))

    #main bridge elements
    for k in 1:nz
        for j = 1:nx
            for i = 1:ny
                bridges = [BridgeElement(primreshaped[i,j,k], x, primreshaped[i+1,j,k], x, joistSection, joistRelease) for x in range(0,1,njoists)[2:end-1]]

                for bridge in bridges
                    bridge.id = :joist
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

    return model;
end