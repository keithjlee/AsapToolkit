function clear_supports!(model::Asap.AbstractModel)

    for node in model.nodes
        fixnode!(node, :free)
        node.id = :free
    end

    process!(model)

end