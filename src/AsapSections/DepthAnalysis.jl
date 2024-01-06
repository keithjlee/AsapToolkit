"""
    intersection(p1::T, p2::T, p3::T, p4::T) where {T <: AbstractVector{Float64}}

Find the intersection point (if it exists) on between two lines:

L1: from p1 to p2
L2: from p3 to p4

returns a boolean and a vector: 
 - (true, point) if an intersection point exists on L1
 - (false, nothing) if an intersection poitn does not exist on L1
"""
function intersection(p1::T, p2::T, p3::T, p4::T) where {T <: AbstractVector{Float64}}

    #expand
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    #intersection parameter on first line (p1 → p2)
    t = ((x1- x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom

    #intersection
    if 0 ≤ t ≤ 1
        return true, p1 + t * (p2 .- p1)
    else
        return false, nothing
    end

end

"""
    sutherland_hodgman(section::PolygonalSection, depth::Float64)

Implements the Sutherland-Hodgman polygon clipping algorithm between a given section and a horizontal line at depth `depth` from the top of the section. "top" is the extreme y-position of the section. Returns a vector of vectors containing the vertices of the new clipped polygon.
"""
function sutherland_hodgman(section::PolygonalSection, depth::Float64; return_section = false)

    #section points
    points = section.points
    npoints = size(points, 2)


    @assert depth ≤ section.ymax - section.ymin

    #clipping edge
    y_clip = section.ymax - depth

    e0 = [section.xmin - 1, y_clip]
    e1 = [section.xmax + 1, y_clip]

    #main algorithm
    new_polygon = Vector{Vector{Float64}}()

    #man algorithm
    for i = 1:npoints
        p1 = section.points_circular[:, i]
        p2 = section.points_circular[:, i + 1]

        if p1[2] ≥ y_clip
            if p2[2] ≥ y_clip
                push!(new_polygon, p2)
            else
                _, point = intersection(p1, p2, e0, e1)
                push!(new_polygon, point)
            end
        else
            if p2[2] ≥ y_clip
                _, point = intersection(p1, p2, e0, e1)
                push!(new_polygon, point)
                push!(new_polygon, p2)
            end
        end
    end

    return_section ? typeof(section)(new_polygon) :  new_polygon

end

"""
    sutherland_hodgman_abs(section::PolygonalSection, y::Float64)

Implements the Sutherland-Hodgman polygon clipping algorithm between a given section and a horizontal line at position `y`. if `y` is above the section, nothing is returned. If `y` is below the section, the points of the original section are returned
"""
function sutherland_hodgman_abs(section::PolygonalSection, y::Float64; return_section = false)

    #section points
    points = section.points
    npoints = size(points, 2)

    # if absolute position is above section, no clip occurs
    if y > section.ymax
        new_polygon = [[0., 0.], [0., 0.]]
    #if absoluteposition encompasses the whole section, return section
    elseif y < section.ymin
        new_polygon = [Vector(col) for col in eachcol(section.points)]
    else
        #clipping edge
        y_clip = y

        e0 = [section.xmin - 1, y_clip]
        e1 = [section.xmax + 1, y_clip]

        #main algorithm
        new_polygon = Vector{Vector{Float64}}()

        #man algorithm
        for i = 1:npoints
            p1 = section.points_circular[:, i]
            p2 = section.points_circular[:, i + 1]

            if p1[2] ≥ y_clip
                if p2[2] ≥ y_clip
                    push!(new_polygon, p2)
                else
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                end
            else
                if p2[2] ≥ y_clip
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                    push!(new_polygon, p2)
                end
            end
        end
    end

    return_section ? typeof(section)(new_polygon) :  new_polygon

end

"""
    depth_map(section::PolygonalSection, n::Integer = 250)

Returns a tuple of (sampled_depths, cumulative_areas) of length n (default = 250) giving the cumulative area covered from the top of the section (maximum y-position) to the bottom of the section (minimum y-position)
"""
function depth_map(section::PolygonalSection, n::Integer = 250)

    #section points
    points = section.points
    npoints = size(points, 2)

    #sampling increments
    depth_range = range(0, section.ymax - section.ymin, n)


    # A = f(depth)

    Astore = Vector{Float64}()

    for depth in depth_range
        y_clip = section.ymax - depth

        e0 = [section.xmin - 1, y_clip]
        e1 = [section.xmax + 1, y_clip]
    
        #main algorithm
        new_polygon = Vector{Vector{Float64}}()
    
        for i = 1:npoints
            p1 = section.points_circular[:, i]
            p2 = section.points_circular[:, i + 1]
    
            if p1[2] ≥ y_clip
                if p2[2] ≥ y_clip
                    push!(new_polygon, p2)
                else
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                end
            else
                if p2[2] ≥ y_clip
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                    push!(new_polygon, p2)
                end
            end
        end

        push!(Astore, poly_area(new_polygon))

    end

    return collect(range(section.ymax, section.ymin, n)), Astore

end

"""
depth_map_abs(section::PolygonalSection, y_positions::Vector{Float64})

Returns a vector of cumulative areas at each point in `y_positions`. If a value `y` in `y_positions` is above the section, a value of 0. is provided. If a value `y` in `y_positions` is below the section, the entire section area is provided.
"""
function depth_map_abs(section::PolygonalSection, y_positions::Vector{Float64})

    #section points
    points = section.points
    npoints = size(points, 2)

    Astore = Vector{Float64}()

    for y in y_positions

        if y > section.ymax
            push!(Astore, 0.)
            continue
        end

        if y < section.ymin
            push!(Astore, section.area)
            continue
        end

        y_clip = y

        e0 = [section.xmin - 1, y_clip]
        e1 = [section.xmax + 1, y_clip]
    
        #main algorithm
        new_polygon = Vector{Vector{Float64}}()
    
        for i = 1:npoints
            p1 = section.points_circular[:, i]
            p2 = section.points_circular[:, i + 1]
    
            if p1[2] ≥ y_clip
                if p2[2] ≥ y_clip
                    push!(new_polygon, p2)
                else
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                end
            else
                if p2[2] ≥ y_clip
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                    push!(new_polygon, p2)
                end
            end
        end

        push!(Astore, poly_area(new_polygon))

    end

    return Astore

end

"""
    depth_map(section::CompoundSection, n::Integer = 250)

Returns the depth map of a compound section
"""
function depth_map(section::CompoundSection, n::Integer = 250)

    #sampling range
    depth_range = collect(range(section.ymax, section.ymin, n))

    #individual depth maps
    depth_maps = [depth_map_abs(sec, depth_range) for sec in section.solids]

    return depth_range, sum(depth_maps)

end


"""
    area_from_depth(section::PolygonalSection, depth::Float64)

Returns the area enclosed by a distance `depth` from the top of a given section
"""
function area_from_depth(section::PolygonalSection, depth::Float64)
    poly_area(sutherland_hodgman(section, depth))
end

"""
    area_from_depth(section::CompoundSection, depth::Float64)

Returns the area enclosed by a distance `depth` from the top of a compound section
"""
function area_from_depth(section::CompoundSection, depth::Float64)
    y = section.ymax - depth

    sum(poly_area.(sutherland_hodgman_abs.(section.solids, y)))
end

"""
    area_from_depth_abs(section::PolygonalSection, y::Float64)

Returns the area enclosed by the intersection of the section and a horizontal line at position `y`. Will return 0 if y is above the section.
"""
function area_from_depth_abs(section::PolygonalSection, y::Float64)
    poly_area(sutherland_hodgman_abs(section, y))
end


"""
    depth_from_area(section::PolygonalSection, area::Float64; max_iter = 500, show_stats = false)

Returns the required depth from the top of the section to achieve a target area.

- `max_iter = 500` sets the maximum number of iterations performed
- `tol = 1e-2` sets the relative stopping tolerance
- `show_stats = true` outputs the error of the solution
"""
function depth_from_area(section::PolygonalSection, area::Float64; max_iter = 500, rel_tol = 1e-3, show_stats = true)

    #first assert target area is possible
    @assert area > 0 && area < section.area "Area must be non-zero and less than the total area of the section"

    #section points
    points = section.points
    npoints = size(points, 2)

    height = section.ymax - section.ymin

    #starting guess at halfway
    depth = height / 2
    upperbound = height
    lowerbound = 0.
    err = 1e3
    iter = 1

    while iter ≤ max_iter
        y_clip = section.ymax - depth

        e0 = [section.xmin - 1, y_clip]
        e1 = [section.xmax + 1, y_clip]
    
        #main algorithm
        new_polygon = Vector{Vector{Float64}}()
    
        for i = 1:npoints
            p1 = section.points_circular[:, i]
            p2 = section.points_circular[:, i + 1]
    
            if p1[2] ≥ y_clip
                if p2[2] ≥ y_clip
                    push!(new_polygon, p2)
                else
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                end
            else
                if p2[2] ≥ y_clip
                    _, point = intersection(p1, p2, e0, e1)
                    push!(new_polygon, point)
                    push!(new_polygon, p2)
                end
            end
        end

        #enclosed area
        solved_area = poly_area(new_polygon)

        #relative error
        err = abs(solved_area - area) / area

        if err < rel_tol
            break
        end

        #difference
        diff = solved_area - area

        if diff < 0
            lowerbound = depth
            depth = (upperbound + depth) / 2
        else
            upperbound = depth
            depth = (lowerbound + depth) / 2
        end
        
        iter += 1

    end

    if show_stats
        println("Iterations: $iter")
        println("Relative Error: $err")
    end

    return depth

end

"""
    depth_from_area(section::PolygonalSection, area::Float64; max_iter = 500, show_stats = false)

Returns the required depth from the top of the section to achieve a target area.

- `max_iter = 500` sets the maximum number of iterations performed
- `tol = 1e-2` sets the relative stopping tolerance
- `show_stats = true` outputs the error of the solution
"""
function depth_from_area(section::CompoundSection, area::Float64; max_iter = 500, rel_tol = 1e-3, show_stats = true)

    #first assert target area is possible
    @assert area > 0 && area < section.area "Area must be non-zero and less than the total area of the section"

    height = section.ymax - section.ymin

    #starting guess at halfway
    depth = height / 2
    upperbound = height
    lowerbound = 0.
    err = 1e3
    iter = 1

    while iter ≤ max_iter
        
        solved_area = area_from_depth(section, depth)

        #difference
        diff = solved_area - area
        
        #relative error
        err = abs(solved_area - area) / area
        if err < rel_tol
            break
        end

        if diff < 0
            lowerbound = depth
            depth = (upperbound + depth) / 2
        else
            upperbound = depth
            depth = (lowerbound + depth) / 2
        end
        
        iter += 1

    end

    if show_stats
        println("Iterations: $iter")
        println("Relative Error: $err")
    end

    return depth

end