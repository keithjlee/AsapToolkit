using Interpolations, kjlMakie, Asap, AsapToolkit, JSON;
set_theme!(kjl_dark)

#sampling parameters
begin
    nxrange = 10:30
    nyrange = 10:30
    dxrange = 1000:250:2000
    dyrange = 1000:250:2000
    dzrange = 1000:250:3500
    supps = [:corner, :x, :y, :xy]
end
# n = 5
# x = range(0, 1, n)
# y = range(0, 1, n)
# z = 3000 .* rand(n,n)

# itp = cubic_spline_interpolation((x,y), z)

# i = range(0,1, 50)
# j = range(0,1, 50)
# k = [itp(i,j) for i in i, j in j]

begin
    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

begin
    pstore = Vector{Vector{Point3{Float64}}}()
    estore = Vector{Vector{Point3{Float64}}}()
    dstore = Vector{Vector{Vector{Float64}}}()
    fstore = Vector{Vector{Float64}}()
    fnstore = Vector{Vector{Float64}}()
    crstore = Vector{Tuple{Float64, Float64}}()
    suppstore = Vector{Vector{Point3{Float64}}}()
end

folder = "spaceframes/"
@time for i = 1:10_000
    println("SAMPLE $i")

    fn = folder * "spaceframe_$i"

    nx = rand(nxrange)
    ny = rand(nyrange)
    dx = rand(dxrange)
    dy = rand(dyrange)
    dz = rand(dzrange)
    supp = rand(supps)

    sf = generatespaceframe(
        nx, 
        dx, 
        ny, 
        dy, 
        dz, 
        tube, 
        ; 
        load = [0., 0., -30e3],
        support = supp
        );

    truss = sf.truss;

    eforces =  getindex.(getproperty.(truss.elements, :forces), 2)
    pts = Point3.(getproperty.(truss.nodes, :position))
    els = vcat([pts[id] for id in getproperty.(truss.elements, :nodeIDs)]...)
    sups = Point3.(getproperty.(truss.nodes[:support], :position))

    if i % 5 == 0
        push!(fstore, eforces)
        push!(fnstore, abs.(eforces) ./ maximum(abs.(eforces)))
        push!(crstore, (-1, 1) .* maximum(abs.(eforces)))
        push!(suppstore, sups)
        push!(estore, els)
        push!(pstore, pts)
        push!(dstore, getproperty.(truss.nodes, :displacement))
    end

    fig = Figure()
    ax = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax); hidespines!(ax)

    linesegments!(els, color = :white)
    scatter!(sups, strokecolor = :white, color = :black)

    save(fn * ".png", fig)


    data = Dict(
        "positions" => [node.position[1:2] for node in truss.nodes],
        "elementindices" => [e.nodeIDs .- 1 for e in truss.elements],
        "supportindices" => findall(truss.nodes, :support) .- 1,
        "compliance" => truss.compliance,
        "loadednodes" => [load.node.nodeID for load in truss.loads] .- 1,
        "loadvalues" => [load.value for load in truss.loads],
        "internalforces" => getindex.(getproperty.(truss.elements, :forces), 2),
        "nodedisplacements" => getproperty.(truss.nodes, :displacement),
        "E" => first(truss.elements).section.E,
        "A" => first(truss.elements).section.A
    )

    ds = JSON.json(data)

    open(fn * ".json","w") do f 
        write(f, ds) 
    end


    fig
end


i = Observable(1)

els = @lift(estore[$i])
sup = @lift(suppstore[$i])

begin
    fig = Figure()
    ax = Axis3(fig[1,1],
        aspect = :data)

    hidedecorations!(ax)
    hidespines!(ax)

    linesegments!(els,
        color = :white)

    scatter!(sup,
        color = :black,
        strokecolor = :white)

    on(i) do _
        reset_limits!(ax)
    end
    fig
end

iterator = 1:20:2000
record(fig, "spaceframes.gif", iterator; framerate = 5) do x
    i[] = x
end