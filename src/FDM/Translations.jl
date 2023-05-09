"""
    toNetwork(model::TrussModel)

Convert a solved truss model into an equivalent FDM Network. 
"""
function toNetwork(model::TrussModel)
    if !model.processed || isnothing(model.u)
        error("Analyze truss model before conversion")
    end

    # convert nodes
    nodeset = Vector{FDMnode}()

    for node in model.nodes
        pos = node.position
        id = node.id

        dof = all(node.dof) ? true : false

        fdmn = FDMnode(pos, dof)
        fdmn.id = id

        push!(nodeset, fdmn)
    end

    #convert loads
    loadset = Vector{FDMload}()

    for load in model.loads
        i = load.node.nodeID

        push!(loadset, nodeset, i, load.value)
    end

    #convert elements
    elset = Vector{FDMelement}()

    for element in model.elements
        istart, iend = element.nodeIDs
        id = element.id
        q = last(element.forces) / element.length

        el = FDMelement(nodeset, istart, iend, q)
        el.id = id

        push!(elset, el)
    end

    network = Network(nodeset, elset, loadset)
    solve!(network)

    return network
end

"""
    toTruss(network::Network, section::AbstractSection)

Convert a solved FDM Network into an equivalent truss model with a given section. All fixed nodes are converted into pinned boundary conditions.
"""
function toTruss(network::Network, section::Asap.AbstractSection)
    if !network.processed
        error("Analyze network before conversion")
    end

    nodeset = Vector{TrussNode}()
    elset = Vector{TrussElement}()
    loadset = Vector{NodeForce}()

    #convert loads
    for node in network.nodes
        pos = [node.x, node.y, node.z]
        id = node.id

        dof = node.dof ? :free : :pinned

        tn = TrussNode(pos, dof)
        tn.id = id

        push!(nodeset, tn)
    end

    #convert elements
    for element in network.elements
        te = TrussElement(nodeset, [element.iStart, element.iEnd], section)

        te.id = element.id

        push!(elset, te)
    end

    #convert loads
    for load in network.loads
        push!(loadset, NodeForce(nodeset[load.index], load.force))
    end

    truss = TrussModel(nodeset, elset, loadset)
    solve!(truss)

    return truss

end