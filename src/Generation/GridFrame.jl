struct GridFrame <: AbstractGenerator
    model::Model
    nx::Integer
    dx::Real
    ny::Integer
    dy::Real
    igrid::Matrix{Int64}

    function GridFrame(nx::Integer, Lx::Real, ny::Integer, Ly::Real, section::Asap.Section; load = [0., 0., -1.], support = :corner, support_type = :pinned)

        dx = Lx / nx
        dy = Ly / ny

        @assert in(support, [:corner, :x, :y, :xy])

        nodal_positions = Vector{Vector{Float64}}()
        index_matrix = zeros(Int64, ny+1, nx+1)

        index = 1
        for i = 1:nx+1
            for j = 1:ny+1
                position = [dx * (j-1), dy * (i-1), 0.]
                push!(nodal_positions, position)

                index_matrix[j, i] = index
                index += 1
            end
        end

        if support == :corner
            support_indices = [index_matrix[1, 1], index_matrix[ny+1, 1], index_matrix[1, nx+1], index_matrix[ny+1, nx+1]]
        elseif support == :x
            support_indices = index_matrix[[1, ny+1], :]
        elseif support == :y
            support_indices = index_matrix[:, [1, nx+1]]
        else
            support_indices = [index_matrix[[1, ny+1], :][:]; index_matrix[2:ny, [1, nx+1]][:]]
        end

        #make nodes
        nodes = [Node(pos, :free, :free) for pos in nodal_positions]

        #make support nodes
        for node in nodes[support_indices]
            fixnode!(node, support_type)
            node.id = :support
        end

        #make elements
        elements = Vector{Element}()

        #horizontal elements
        for i in axes(index_matrix, 1)
            for j = 1:nx
                
                node_indices = index_matrix[i, [j, j+1]]

                push!(elements, Element(nodes, node_indices, section))
            end
        end

        #vertical elements
        for j in axes(index_matrix, 2)
            for i = 1:ny
                node_indices = index_matrix[[i, i+1], j]

                push!(elements, Element(nodes, node_indices, section))
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
            index_matrix
        )
    end

end