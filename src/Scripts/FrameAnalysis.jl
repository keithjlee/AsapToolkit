using Asap, AsapToolkit, kjlMakie
set_theme!(kjl_dark)

# choose sections
begin
    W = allW()
    iStiffness = sortperm(getproperty.(W, :Ix))
    nsecs = length(W)

    mat = Steel_Nmm
    girder = toASAPframe((W[iStiffness])[Int(round(nsecs * 0.5))], mat.E, mat.G)
    joist = toASAPframe((W[iStiffness])[Int(round(nsecs * .2))], mat.E, mat.G)

    iArea = sortperm(getproperty.(W, :A))
    col = toASAPframe((W[iArea])[Int(round(nsecs * .75))], mat.E, mat.G)

    H = allHSSRound()
    iArea = sortperm(getproperty.(H, :A))
    nh = length(H)
    tube = toASAPframe((H[iArea])[Int(round(nh * 0.9))], mat.E, mat.G)

    girder.ρ = joist.ρ = tube.ρ = mat.ρ
end

# geometry hyperparameters for frame
begin
    nx = 3
    ny = 6
    nz = 30
    dx = 6000
    dy = 5000
    dz = 4200
end

#analysis
begin
    @time frame = AsapToolkit.generateframe(nx,
        dx,
        ny,
        dy,
        nz,
        dz,
        500,
        col,
        girder,
        joist,
        tube;
        columnPsi = pi/2,
        diaphragm = true,
        );

    model = frame.model;
    loads = [LineLoad(j, [0., 0., -30.]) for j in model.elements[:joist]];
    @time solve!(model, loads)

    dispfac = Observable(2.)
    resolution = 20

    ixn = frame.iExteriorXnodes
    iyn = frame.iExteriorYnodes
    ixj = frame.iExteriorXjoists
    iyp = frame.iExteriorYprimaries
end

#optional new loads
begin
    windloads = [NodeForce(n, [0., 40e3, 0.]) for n in model.nodes[ixn]]
    @time solve!(model, [loads; windloads])
end

#meta analysis
increment = 50.
begin
    forceresults = forces(model, increment)
    dispresults = displacements(model, increment)


    e_ids = [fr.element.id for fr in forceresults]
    ijoist = findall(e_ids .== :joist)
    iprim = findall(e_ids .== :primary)
    icol = findall(e_ids .== :column)

    unitDict = Dict(:P => "[N]",
        :Vy => "[N]",
        :Vz => "[N]",
        :My => "[Nmm]",
        :Mz => "[Nmm]")

end

#plotting elements
begin
    p_u = Point3.(getproperty.(model.nodes, :position))
    e_u = vcat([p_u[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    p_d = @lift(p_u .+ $dispfac .* [n.displacement[1:3] for n in model.nodes])
    e_d = [@lift(Point3.(eachcol(d.basepositions .+ $dispfac .* d.uglobal))) for d in dispresults]

    # element wise values
    crfactor = Observable(0.5)
    lwfactor = Observable(5.0)
    property = Observable(:P)

    propertyvals = [@lift(getproperty(fr, $property)) for fr in forceresults]
    propmax = @lift(maximum([maximum(abs.(getproperty(fr, $property))) for fr in forceresults]))
    propertyvalsnormalized = [@lift($lwfactor .* abs.(getproperty(fr, $property)) ./ $propmax .+ 2) for fr in forceresults]
    propertyrange = @lift($crfactor .* (-1, 1) .* $propmax)
    proprangefull = @lift((-1, 1) .* $propmax)
    propname = @lift(string($property) * " " * unitDict[$property])
end

# structure plot
begin
    fig = Figure(resolution = (1000,1000))
    ax = Axis3(fig[1,1],
        protrusions = 100,
        zlabeloffset = 75,

        aspect = :data)

    gridtoggle!(ax); simplifyspines!(ax)

    eu = linesegments!(e_u,
        color = :white,
        linewidth = .5,
        )

    eu.visible = false

    ed = [lines!(ed,
        color = prop,
        colorrange = propertyrange,
        colormap = pink2blue,
        linewidth = propvaln
        ) for (ed, prop, propvaln) in zip(e_d, propertyvals, propertyvalsnormalized)]

    cb = Colorbar(fig[1,2],
        # vertical = false,
        # flipaxis = false,
        colorrange = proprangefull,
        colormap = pink2blue,
        label = propname)

    fig
end

#force analysis
begin
    xstore = [x .- x[end] / 2 for x in getproperty.(forceresults, :x)]
    my = [my ./ 1e6 for my in getproperty.(forceresults, :My)]
    vy = [vy ./ 1e3 for vy in getproperty.(forceresults, :Vy)]
    mz = [mz ./ 1e6 for mz in getproperty.(forceresults, :Mz)]
    vz = [vz ./ 1e3 for vz in getproperty.(forceresults, :Vz)]
    p = [p ./ 1e3 for p in getproperty.(forceresults, :P)]

    uids = unique(e_ids)
    nuids = length(uids)
    coldict = Dict(uids .=> collect(1:nuids))

    cm = :Wistia;
    cg = discretize(nuids + 2, colormap = cm)[2:end-1];

    colnums = [cg[coldict[id]] for id in e_ids]

    legelems = [LineElement(color = col) for col in [cg[coldict[id]] for id in uids]]
    lw = 3
end

#visualize force distribution
begin
    forcefig = Figure(backgroundcolor = :black)

    ax_my = Axis(forcefig[1,1],
        aspect = nothing,
        ylabel = "My [kNm]")
    
    hidexdecorations!(ax_my)

    hlines!(ax_my, [0.], color = :white)
    M = [lines!(x, m,
        color = (col, 0.5),
        linewidth = lw,
        ) for (x, m, col) in zip(xstore, my, colnums)]
    
    ax_vy = Axis(forcefig[2,1],
        aspect = nothing,
        xlabel = "x [mm]",
        ylabel = "V [kN]")

    hlines!(ax_vy, [0.], color = :white)
    V =[lines!(x, v,
        color = (col, 0.5),
        linewidth = lw
        ) for (x, v, col) in zip(xstore, vy, colnums)]

    ax_p = Axis(forcefig[1:2,2],
        aspect = nothing,
        yaxisposition = :right,
        xlabel = "P [kN]",
        ylabel = "x [mm]")

    vlines!(ax_p, [0.], color = :white)
    P = [lines!(pp, x,
        color = (col, 0.5),
        linewidth = lw
        ) for (x,pp,col) in zip(xstore, p, colnums)]


    leg = Legend(forcefig[1:2, 3],
        legelems,
        string.(uids))

    forcefig
end

#as scattered points
begin
    forcefig = Figure(backgroundcolor = :black)

    ax_my = Axis(forcefig[1,1],
        aspect = nothing,
        ylabel = "My [kNm]")
    
    hidexdecorations!(ax_my)

    hlines!(ax_my, [0.], color = :white)
    M = [scatter!(x, m,
        color = (col, 0.5),
        # linewidth = lw,
        ) for (x, m, col) in zip(xstore, my, colnums)]
    
    ax_vy = Axis(forcefig[2,1],
        aspect = nothing,
        xlabel = "x [mm]",
        ylabel = "V [kN]")

    hlines!(ax_vy, [0.], color = :white)
    V =[scatter!(x, v,
        color = (col, 0.5),
        # linewidth = lw
        ) for (x, v, col) in zip(xstore, vy, colnums)]

    ax_p = Axis(forcefig[1:2,2],
        aspect = nothing,
        yaxisposition = :right,
        xlabel = "P [kN]",
        ylabel = "x [mm]")

    vlines!(ax_p, [0.], color = :white)
    P = [scatter!(pp, x,
        color = (col, 0.5),
        # linewidth = lw
        ) for (x,pp,col) in zip(xstore, p, colnums)]


    leg = Legend(forcefig[1:2, 3],
        legelems,
        string.(uids))

    forcefig
end