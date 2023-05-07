export MfuncDict
# MfuncDict = Dict((LineLoad, :fixedfixed) => MLine_fixedfixed,
#     (LineLoad, :freefree) => MLine_freefree,
#     (LineLoad, :fixedfree) => MLine_fixedfree,
#     (LineLoad, :freefixed) => MLine_freefixed,
#     (LineLoad, :joist) => MLine_freefree,
#     (PointLoad, :fixedfixed) => MPoint_fixedfixed,
#     (PointLoad, :freefree) => MPoint_freefree,
#     (PointLoad, :fixedfree) => MPoint_fixedfree,
#     (PointLoad, :freefixed) => MPoint_freefixed,
#     (PointLoad, :joist) => MPoint_freefree
#     )

MLineLoad = Dict(:fixedfixed => MLine_fixedfixed,
    :freefree => MLine_freefree,
    :fixedfree => MLine_fixedfree,
    :freefixed => MLine_freefixed,
    :joist => MLine_freefree
    )

MPointLoad = Dict(:fixedfixed => MPoint_fixedfixed,
    :freefree => MPoint_freefree,
    :fixedfree => MPoint_fixedfree,
    :freefixed => MPoint_freefixed,
    :joist => MPoint_freefree)

export VfuncDict
# VfuncDict = Dict((LineLoad, :fixedfixed) => VLine_fixedfixed,
#     (LineLoad, :freefree) => VLine_freefree,
#     (LineLoad, :fixedfree) => VLine_fixedfree,
#     (LineLoad, :freefixed) => VLine_freefixed,
#     (LineLoad, :joist) => VLine_freefree,
#     (PointLoad, :fixedfixed) => VPoint_fixedfixed,
#     (PointLoad, :freefree) => VPoint_freefree,
#     (PointLoad, :fixedfree) => VPoint_fixedfree,
#     (PointLoad, :freefixed) => VPoint_freefixed,
#     (PointLoad, :joist) => VPoint_freefree)

VLineLoad = Dict(:fixedfixed => VLine_fixedfixed,
    :freefree => VLine_freefree,
    :fixedfree => VLine_fixedfree,
    :freefixed => VLine_freefixed,
    :joist => VLine_freefree)

VPointLoad = Dict(:fixedfixed => VPoint_fixedfixed,
    :freefree => VPoint_freefree,
    :fixedfree => VPoint_fixedfree,
    :freefixed => VPoint_freefixed,
    :joist => VPoint_freefree)

export DfuncDict
# DfuncDict = Dict((LineLoad, :fixedfixed) => DLine_fixedfixed,
#     (LineLoad, :freefree) => DLine_freefree,
#     (LineLoad, :fixedfree) => DLine_fixedfree,
#     (LineLoad, :freefixed) => DLine_freefixed,
#     (LineLoad, :joist) => DLine_freefree,
#     (PointLoad, :fixedfixed) => DPoint_fixedfixed,
#     (PointLoad, :freefree) => DPoint_freefree,
#     (PointLoad, :fixedfree) => DPoint_fixedfree,
#     (PointLoad, :freefixed) => DPoint_freefixed,
#     (PointLoad, :joist) => DPoint_freefree)

DLineLoad = Dict(:fixedfixed => DLine_fixedfixed,
    :freefree => DLine_freefree,
    :fixedfree => DLine_fixedfree,
    :freefixed => DLine_freefixed,
    :joist => DLine_freefree)

DPointLoad = Dict(:fixedfixed => DPoint_fixedfixed,
    :freefree => DPoint_freefree,
    :fixedfree => DPoint_fixedfree,
    :freefixed => DPoint_freefixed,
    :joist => DPoint_freefree)
    
export release2DOF
release2DOF = Dict(:fixedfixed => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    :freefixed => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    :fixedfree => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    :freefree => [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    :joist => [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0])

export planarDOFs
planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])



function accumulate!(load::LineLoad, 
    xvals::Vector{Float64}, 
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    release = load.element.release

    # distributed load magnitudes in LCS
    wy, wz = - (R  * load.value)[2:3]

    # extract relevant function
    mfunction = MLineLoad[release]
    vfunction = VLineLoad[release]

    My .+= mfunction.(wy, L, xvals)
    Vy .+= vfunction.(wy, L, xvals)

    Mz .+= mfunction.(wz, L, xvals)
    Vz .+= vfunction.(wz, L, xvals)
end

function accumulate!(load::PointLoad, 
    xvals::Vector{Float64}, 
    My::Vector{Float64}, 
    Vy::Vector{Float64}, 
    Mz::Vector{Float64}, 
    Vz::Vector{Float64})

    R = load.element.R[1:3, 1:3]
    L = load.element.length
    release = load.element.release
    frac = load.position

    # distributed load magnitudes in LCS
    py, pz = - (R  * load.value)[2:3]

    # extract relevant function
    mfunction = MPointLoad[release]
    vfunction = VPointLoad[release]

    My .+= mfunction.(py, L, xvals, frac)
    Vy .+= vfunction.(py, L, xvals, frac)

    Mz .+= mfunction.(pz, L, xvals, frac)
    Vz .+= vfunction.(pz, L, xvals, frac)
end

struct ElementResults
    element::Element
    resolution::Integer
    x::Vector{Float64}
    My::Vector{Float64}
    Vy::Vector{Float64}
    Mz::Vector{Float64}
    Vz::Vector{Float64}
end

export internalforces
function internalforces(element::Element, model::Model; resolution = 20)
    
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
    Vystart, Mystart, Vzstart, Mzstart = Flocal[[2, 6, 3, 5]] .* [1, 1, 1, -1]

    # initialize internal force vectors
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in model.loads[element.loadIDs]
        accumulate!(load,
            xinc,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return ElementResults(element, resolution, xinc, My, Vy, Mz, Vz)
end

function internalforces(element::Element, loads::Vector{<:Asap.ElementLoad}; resolution = 20)
    
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
    Vystart, Mystart, Vzstart, Mzstart = Flocal[[2, 6, 3, 5]] .* [1, 1, 1, -1]

    # initialize internal force vectors
    My = Vystart .* xinc .- Mystart
    Vy = zero(My) .+ Vystart
    Mz = Vzstart .* xinc .- Mzstart
    Vz = zero(Mz) .+ Vzstart

    # accumulate loads
    for load in loads[element.loadIDs]
        accumulate!(load,
            xinc,
            My,
            Vy,
            Mz,
            Vz)  
    end

    return ElementResults(element, resolution, xinc, My, Vy, Mz, Vz)
end

function internalforces(elements::Vector{<:Asap.FrameElement}, model::Model; resolution = 20)
    
    xstore = Vector{Float64}()
    mystore = Vector{Float64}()
    vystore = Vector{Float64}()
    mzstore = Vector{Float64}()
    vzstore = Vector{Float64}()

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
        Vystart, Mystart, Vzstart, Mzstart = Flocal[[2, 6, 3, 5]] .* [1, 1, 1, -1]

        # initialize internal force vectors
        My = Vystart .* xinc .- Mystart
        Vy = zero(My) .+ Vystart
        Mz = Vzstart .* xinc .- Mzstart
        Vz = zero(Mz) .+ Vzstart

        # accumulate loads
        for load in model.loads[element.loadIDs]
            accumulate!(load,
                xinc,
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

        mystore = [mystore; My]
        vystore = [vystore; Vy]
        mzstore = [mzstore; Mz]
        vzstore = [vzstore; Vz]

    end

    return ElementResults(elements[1], resolution, xstore, mystore, vystore, mzstore, vzstore)
end
