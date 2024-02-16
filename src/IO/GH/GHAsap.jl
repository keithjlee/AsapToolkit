include("GHNode.jl")
include("GHElements.jl")
include("GHLoad.jl")
include("GHModel.jl")

"""
    GHsave(model::Asap.AbstractModel, filename::String)

Save an AsapModel as a .json file with GHmodel/GHnode/GHelement/GHload data structure. Add ".json" to the end of `filename` is optional.
"""
function GHsave(model::Asap.AbstractModel, filename::String)
    if filename[end-4:end] != ".json"
        filename *= ".json"
    end

    ghmodel = GHmodel(model)
    open(filename, "w") do f
        write(f, JSON.json(ghmodel))
    end
end

"""
    GHsave(model::GHmodel, filename::String)

Save a GHmodel as a .json file with GHmodel/GHnode/GHelement/GHload data structure. Add ".json" to the end of `filename` is optional.
"""
function GHsave(model::GHmodel, filename::String)
    if filename[end-4:end] != ".json"
        filename *= ".json"
    end

    open(filename, "w") do f
        write(f, JSON.json(model))
    end
end