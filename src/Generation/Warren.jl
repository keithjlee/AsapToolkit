struct Warren2D <: AbstractGenerator
    model::TrussModel
    n::Integer
    dx::Real
    dy::Real
    section::Asap.AbstractSection
    type::Symbol
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
"""
function Warren2D(n::Integer,
        dx::Real,
        dy::Real,
        section::Asap.AbstractSection;
        load = [0., -1., 0.],
        type = :arch)

    @assert n % 2 != 0 "n must be odd"
    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{TrussNode}()

    #generate longer chord first
    if type == :arch
        longid = :bottomchord
        shortid = :topchord
        y = dy
    else
        longid = :topchord
        shortid = :bottomchord
        y = -dy
    end

    for i = 1:n
        xposition = dx * (i - 1)

        node = TrussNode([xposition, 0., 0.], :free)
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

        node = TrussNode([xposition, y, 0.], :free)
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
    loads = [NodeForce(n, load) for n in nodes[longid]]

    #assemble and solve
    model = TrussModel(nodes, elements, loads)
    planarize!(model)
    solve!(model)

    #collect data
    truss = Warren2D(model, n, dx, dy, section, type)

    #output
    return truss
end

function Warren2D(xpositions::Vector{<:Real},
    ypositions::Vector{<:Real},
    ypositions2::Vector{<:Real},
    section::Asap.AbstractSection;
    type = :arch)

    @assert length(xpositions) == length(ypositions) == length(ypositions2) + 1
    @assert type == :arch || type == :catenary "type must be :arch or :catenary"

    #counters
    count = 1
    longids = Vector{Int64}()

    #node collector
    nodes = Vector{TrussNode}()

    #generate longer chord first
    if type == :arch
        longid = :bottomchord
        shortid = :topchord
    else
        longid = :topchord
        shortid = :bottomchord
    end

    n = length(xpositions)
    i = 1

    ## generate long chord up to symmetry
    for (x, y) in zip(xpositions, ypositions)

        node = TrussNode([x, y, 0.], :free)
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
    incs = Lhalf .- reverse(xpositions[1:end-1])
    for (inc, y) in zip(incs, reverse(ypositions[1:end-1]))
        node = TrussNode([inc, y, 0.], :free)
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

        node = TrussNode([xposition, ypositions2[i], 0.], :free)
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
    truss = Warren2D(model, n, dx, dy, section, type)

end