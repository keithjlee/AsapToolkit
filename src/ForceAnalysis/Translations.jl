const etype2DOF = Dict(
    Element{Asap.FixedFixed} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    Element{Asap.FreeFree} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1],
    Element{Asap.FixedFree} => [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
    Element{Asap.FreeFixed} => [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
    Element{Asap.Joist} => [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0]
)

export planarDOFs
const planarDOFs = Dict(:X => [1, 7],
    :XY => [2, 6, 8, 12],
    :XZ => [3, 5, 9, 11])

"""
Collect all elements that are part of the same continuous member (ie collect shattered elements)
"""
function groupbyid(elements::Vector{<:Asap.AbstractElement})
    ichecked = Vector{Int64}()
    indices = Vector{Vector{Int64}}()

    elementids = getproperty.(elements, :elementID)

    for id in elementids
        
        in(id, ichecked) && continue

        igroup = findall(elementids .== id)

        push!(indices, igroup)
        push!(ichecked, id)
    end

    return indices
end