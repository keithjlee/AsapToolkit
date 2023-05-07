"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulateforce!(load::LineLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    release = load.element.release

    # distributed load magnitudes in LCS
    wx, wy, wz = (R  * load.value) .* [1, -1, -1]

    # extract relevant function
    mfunction = MLineLoad[release]
    vfunction = VLineLoad[release]

    P .+= PLine.(wx, L, xvals)

    My .+= mfunction.(wy, L, xvals)
    Vy .+= vfunction.(wy, L, xvals)

    Mz .+= mfunction.(wz, L, xvals)
    Vz .+= vfunction.(wz, L, xvals)
end

"""
Accumlate the internal forces cause by a given load to the current element
"""
function accumulateforce!(load::PointLoad, 
    xvals::Vector{Float64}, 
    P::Vector{Float64},
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    release = load.element.release
    frac = load.position

    # distributed load magnitudes in LCS
    px, py, pz = (R  * load.value) .* [1, -1, -1]

    # extract relevant function
    mfunction = MPointLoad[release]
    vfunction = VPointLoad[release]

    P .+= PPoint.(px, L, xvals, frac)

    My .+= mfunction.(py, L, xvals, frac)
    Vy .+= vfunction.(py, L, xvals, frac)

    Mz .+= mfunction.(pz, L, xvals, frac)
    Vz .+= vfunction.(pz, L, xvals, frac)
end

"""
Store the internal force results for a given element
"""
struct InternalForces
    element::Element
    resolution::Integer
    x::Vector{Float64}
    P::Vector{Float64}
    My::Vector{Float64}
    Vy::Vector{Float64}
    Mz::Vector{Float64}
    Vz::Vector{Float64}
end


"""
Internal force sampling for an element 
"""
function InternalForces(element::Element, model::Model; resolution = 20)
    
    #beam information
    release = element.release
    L = element.length

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* release2DOF[release]

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in model.loads[element.loadIDs]
        accumulateforce!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return InternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end

"""
    InternalForces(element::Element, loads::Vector{<:ElementLoad}; resolution = 20)

Get internal force results for a given element from a set of loads
"""
function InternalForces(element::Element, loads::Vector{<:Asap.ElementLoad}; resolution = 20)
    
    #beam information
    release = element.release
    L = element.length

    #discretization
    xinc = collect(range(0, L, resolution))

    #end node information
    uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]
    
    # end forces that are relevant to the given element/release condition
    Flocal = (element.R * element.K * uglobal) .* release2DOF[release]

    # shear/moment acting at the *starting* point of an element in LCS
    Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

    # initialize internal force vectors
    P = repeat([Pstart], resolution)
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in loads[element.loadIDs]
        accumulateforce!(load,
            xinc,
            P,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return InternalForces(element, resolution, xinc, P, My, Vy, Mz, Vz)
end

"""
    InternalForces(element::Vector{<:FrameElement}, model::Model; resolution = 20)

Get internal force results for a group of ordered elements that form a single physical element
"""
function InternalForces(elements::Vector{<:Asap.FrameElement}, model::Model; resolution = 20)
    
    xstore = Vector{Float64}()
    pstore = Vector{Float64}()
    mystore = Vector{Float64}()
    vystore = Vector{Float64}()
    mzstore = Vector{Float64}()
    vzstore = Vector{Float64}()

    resolution = Int(round(resolution / length(elements)))

    #beam information
    for element in elements
        release = element.release
        L = element.length

        #discretization
        xinc = collect(range(0, L, resolution))

        #end node information
        uglobal = [element.nodeStart.displacement; element.nodeEnd.displacement]
        
        # end forces that are relevant to the given element/release condition
        Flocal = (element.R * element.K * uglobal) .* release2DOF[release]

        # shear/moment acting at the *starting* point of an element in LCS
        Pstart, Vystart, Mystart, Vzstart, Mzstart = Flocal[[1, 2, 6, 3, 5]] .* [-1, 1, 1, 1, -1]

        # initialize internal force vectors
        P = repeat([Pstart], resolution)
        My = Vystart .* xinc .- Mystart
        Vy = zero(My) .+ Vystart
        Mz = Vzstart .* xinc .- Mzstart
        Vz = zero(Mz) .+ Vzstart

        # accumulate loads
        for load in model.loads[element.loadIDs]
            accumulateforce!(load,
                xinc,
                P,
                My,
                Vy,
                Mz,
                Vz)  
        end

        if isempty(xstore)
            xstore = [xstore; xinc]
        else
            xstore = [xstore; xstore[end] .+ xinc]
        end

        pstore = [pstore; P]
        mystore = [mystore; My]
        vystore = [vystore; Vy]
        mzstore = [mzstore; Mz]
        vzstore = [vzstore; Vz]

    end

    return InternalForces(elements[1], resolution, xstore, pstore, mystore, vystore, mzstore, vzstore)
end



"""
    forceanalysis(model::Model, increment::Real)

Get a vector of element internal force results every `increment`
"""
function forceanalysis(model::Model, increment::Real)

    results = Vector{InternalForces}()

    ids = groupbyid(model.elements)

    for id in ids
        elements = model.elements[id]
        L = sum(getproperty.(elements, :length))
        n = max(Int(round(L / increment)), 2)

        push!(results, InternalForces(elements, model; resolution = n))
    end

    return results
end