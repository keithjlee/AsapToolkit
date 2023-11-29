struct GridFrame <: AbstractGenerator
    model::Model
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    igrid::Matrix{Int64}

    function GridFrame(Lx::Real, nx::Integer, Ly::Real, ny::Integer, section::Asap.Section; load = [0., 0., -1.], support = :corner, support_type = :pinned)

        @assert in(support, [:corner, :x, :y, :xy])

        x_positions = range(0, Lx, nx)
        y_positions = range(0, Ly, ny)

        dx = Lx / (nx-1)
        dy = Ly / (ny-1)

        xyz = Vector{Vector{Float64}}()
        Xmatrix = zeros(ny, nx)
        Ymatrix = zeros(ny, nx)
        igrid = zeros(Int64, ny, nx)

        index = 1
        for iy = 1:nx
            for ix = 1:ny

                x = x_positions[iy]
                y = y_positions[ix]

                igrid[ix, iy] = index
                index += 1

                push!(xyz, [x, y, 0.])
                Xmatrix[ix, iy] = x
                Ymatrix[ix, iy] = y

            end
        end

        if support == :corner
            support_indices = [igrid[1, 1], igrid[ny, 1], igrid[1, nx], igrid[ny, nx]]
        elseif support == :x
            support_indices = igrid[[1, ny], :]
        elseif support == :y
            support_indices = igrid[:, [1, nx]]
        else
            support_indices = [igrid[[1, ny], :][:]; igrid[2:ny-1, [1, nx]][:]]
        end

        #make nodes
        nodes = [Node(pos, :free, :free) for pos in xyz]

        #make support nodes
        for node in nodes[support_indices]
            fixnode!(node, support_type)
            node.id = :support
        end

        #make elements
        elements = Vector{Element}()

        #horizontal elements
        for i = 1:ny
            for j = 1:nx-1
                index = [igrid[i,j], igrid[i,j+1]]
                push!(elements, Element(nodes, index, section))
            end
        end

        #vertical elements
        for j = 1:nx
            for i = 1:ny-1
                index = [igrid[i,j], igrid[i+1,j]]
                push!(elements, Element(nodes, index, section)) 
            end
        end

        #loads
        loads = [NodeForce(node, load) for node in nodes[:free]]

        #assemble
        model = Model(nodes, elements, loads)
        solve!(model)

        new(
            model,
            nx,
            dx,
            ny,
            dy,
            igrid
        )
    end

end