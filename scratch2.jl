begin
    n1 = Node([0., 0., 0.], :pinned)
    n1b = Node([4000., 0., 0.], :free)
    n2 = Node([8000., 0., 0.], :fixed)

    p1 = Element(n1, n1b, girder)
    p2 = Element(n1b, n2, girder)
    # p.Ψ = 0.

    l1 = LineLoad(p1, [0., 0., -20.])
    l2 = LineLoad(p2, [0., 0., -20.])
    l3 = PointLoad(p1, .1, [0., 0., -30e3])
    l4 = PointLoad(p2, 0.5, [0., 0., 50e3])
    l5 = PointLoad(p2, 0.7, [0., 0., -75e3])

    nodes = [n1, n1b, n2]; elements = [p1, p2]; loads = [l1, l2, l3, l4, l5]

    model = Model(nodes, elements, loads)
    solve!(model)

    df = Observable(1000.)

    xyz1 = @lift(Point3.(eachcol(displacedshape(p1; n = n, factor = $df))))
    xyz2 = @lift(Point3.(eachcol(displacedshape(p2; n = n, factor = $df))))

    undisp = Point3.([endpoints(p1); endpoints(p2)])
end

begin
    fig = Figure(backgroundcolor = :black)
    ax = Axis3(fig[1,1],
        aspect = :data)

    simplifyspines!(ax)

    ud = lines!(undisp,
        color = :white,
        linestyle = :dash)

    lines!(xyz1, color = :white, linewidth = 3)
    lines!(xyz2, color = :white, linewidth = 3)

    on(df) do _ reset_limits!(ax) end

    fig
end

resolution = 100

disp_dof_LCS = Asap.localdisplacements(p2; n= resolution)
uglobal = [p1.nodeStart.displacement; p1.nodeEnd.displacement]
force_dof_LCS = p1.R * (p1.K * uglobal + p.Q)
Fy = force_dof_LCS[2]
My = force_dof_LCS[6]

V_dof = Fy
M_dof = Fy * xincs .- My

# DISPLACEMENT FROM LOAD
loads_on_e = model.loads[p2.loadIDs]

load = loads_on_e[1]

Mlocaly = zeros(resolution)
Vlocaly = zeros(resolution)
Dlocaly = zeros(resolution)

for load in loads_on_e
    plocal = p.R[1:3,1:3] * load.value
    w_localy = -plocal[2]
    w_localz = -plocal[3]

    xincs = range(0, p.length, resolution)

    p.release = :fixedfixed
    mfunction = MfuncDict[typeof(load), p.release]
    try
        Mlocaly += mfunction.(w_localy, p.length, xincs)
    catch
        Mlocaly += mfunction.(w_localy, p.length, xincs, load.position)
    end

    vfunction = VfuncDict[typeof(load), p.release]
    
    try
        Vlocaly += vfunction.(w_localy, p.length, xincs)
    catch
        Vlocaly += vfunction.(w_localy, p.length, xincs, load.position)
    end

    dfunction = DfuncDict[typeof(load), p.release]

    try
        Dlocaly += -dfunction.(w_localy, p.length, xincs, p.section.E, p.section.Izz)
    catch
        Dlocaly += -dfunction.(w_localy, p.length, xincs, load.position, p.section.E, p.section.Izz)

    end
end

begin
    fig = Figure(backgroundcolor = :black)

    axM = Axis(fig[1,1], yreversed = true, aspect = nothing)

    lines!(xincs, Mlocaly)

    axV = Axis(fig[2,1], aspect = nothing)

    lines!(xincs, Vlocaly)

    axD = Axis(fig[3,1], aspect = nothing)

    lines!(xincs, Dlocaly)


    hlines!.((axM, axV, axD), [0.], color = :white)

    fig
end


unodal = model.S[model.freeDOFs, model.freeDOFs] \ model.P[model.freeDOFs]
