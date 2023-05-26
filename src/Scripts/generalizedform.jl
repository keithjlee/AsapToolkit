using Asap, AsapToolkit
using Zygote, LinearAlgebra, FiniteDiff
using kjlMakie; set_theme!(kjl_dark)

### Create a spaceframe

#meta parameters
begin
    nx = 25
    dx = 750.
    ny = 15
    dy = 1000.
    dz = 1500.

    sec = rand(allHSSRound())
    tube = toASAPtruss(sec, Steel_Nmm.E)
    σ = 350.
end

n = 5
x = range(0, 1, n)
y = range(0, 1, n)
z = rand(n,n) .* 3000

using Interpolations
itp = cubic_spline_interpolation((x,y), z)

#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, itp, tube, true; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

begin

    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in getproperty.(model.elements, :nodeIDs)]...)

    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    

    linesegments!(e0, color = :white)


    fig
end

#make variables
vars = Vector{AsapToolkit.TrussVariables}()

# nodal Z variables
lb = -1000.
ub = 3500.

lbb = -3500.
ubb = 1000.

#nodal xy variables for bottom nodes
lbxy = -750.
ubxy = 750.

for node in model.nodes
    if node.id == :top
        push!(vars, SpatialVariable(node, 0., lb, ub, :Z))
    end

    if node.id == :bottom
        push!(vars, SpatialVariable(node, 0., lbb, ubb, :Z))
    end

    if node.id == :bottom
        push!(vars, SpatialVariable(node, 0., lbxy, ubxy, :X))
        push!(vars, SpatialVariable(node, 0.,  lbxy, ubxy, :Y))
    end
end

# All bottom elements at once
bottomAreaMaster = AreaVariable(model.elements[sf.ibottom[1]], 500., 0., 35000.)
push!(vars, bottomAreaMaster)

for i in sf.ibottom[2:end]
    push!(vars, CoupledVariable(model.elements[i], bottomAreaMaster))
end

# individual area optimization for web
for element in model.elements[:web]
    push!(vars, AreaVariable(element, 500., 0., 20000.))
end

@time problem = TrussOptProblem(model, vars);

function obj(values::Vector{Float64}, p::TrussOptProblem)
    indexer = p.indexer

    Xnew = addvalues(p.X, indexer.iX, values[indexer.iXg])
    Ynew = addvalues(p.Y, indexer.iY, values[indexer.iYg])
    Znew = addvalues(p.Z, indexer.iZ, values[indexer.iZg])

    Anew = replacevalues(p.A, indexer.iA, values[indexer.iAg])

    ks = [kglobal(Xnew, Ynew, Znew, e, a, id) for (e, a, id) in zip(p.E, Anew, p.params.nodeids)]

    K = assembleglobalK(ks, p.params)

    U = solveU(K, p.params)

    U' * p.params.P[p.params.freeids]
end

vals = problem.values
@time Zygote.gradient(var -> obj(var, problem), vals)[1];