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

#generation and extraction
sf = generatespaceframe(nx, dx, ny, dy, dz, tube; load = [0., 0., -30e3], support = :xy);
model = sf.truss;

# random new loads
newloads = Vector{NodeForce}()
for _ = 1:10
    iset = rand(sf.isquares)
    for i in iset
        push!(newloads, NodeForce(model.nodes[i], [0., 0., -35e3]))
    end
end

solve!(model, [model.loads; newloads])


struct TestParams <: AsapToolkit.AbstractParams
    nodeids::Vector{Vector{Int64}} # [[iNodeStart, iNodeEnd] for element in elements]
    dofids::Vector{Vector{Int64}} # [[dofStartNode..., dofEndNode...] for element in elements]
    P::Vector{Float64} # External load vector
    freeids::Vector{Int64} # [DofFree1, DofFree2,...]
    inzs::Vector{Vector{Int64}} # Indices of elemental K in global S.nzval
    n::Int64 #number of DOF total
    ndofe::Int64 #number of DOF in elemental stiffness matrix
    cp::Vector{Int64} #S.colptr
    rv::Vector{Int64} #S.rowval
    nnz::Int64 #length(S.nzval)
    E::Vector{Float64}
    A::Vector{Float64}
    X::Vector{Float64}
    Y::Vector{Float64}
    Z::Vector{Float64}
end

function TestParams(model::TrussModel)
    nodeids = getproperty.(model.elements, :nodeIDs)
    dofids = getproperty.(model.elements, :globalID)
    P = model.P
    freeids = model.freeDOFs
    inzs = allinz(model)
    n = model.nDOFs
    ndofe = 6
    cp = model.S.colptr
    rv = model.S.rowval
    nnz = length(model.S.nzval)
    E = getproperty.(getproperty.(model.elements, :section), :E)
    A = getproperty.(getproperty.(model.elements, :section), :A)
    X = getindex.(getproperty.(model.nodes, :position), 1)
    Y = getindex.(getproperty.(model.nodes, :position), 2)
    Z = getindex.(getproperty.(model.nodes, :position), 3)

    return TestParams(nodeids, dofids, P, freeids, inzs, n, ndofe, cp, rv, nnz, E, A, X, Y, Z)
end

p = TestParams(model)

vars = Vector{AsapToolkit.OptVariable}()

# nodal Z variables
lb = -1000.
ub = 3500.

lbb = -3500.
ubb = 1000.

#nodal xy variables for bottom nodes
lbxy = -500.
ubxy = 500.

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

# area of bottom elements
for element in model.elements
    if element.id == :bottom
        push!(vars, InternalVariable(element, 500., 0., 35000.))
    end
end

# extract collected values
lowerbounds = getproperty.(vars, :lb)
upperbounds = getproperty.(vars, :ub)
initialvalues = getproperty.(vars, :val)

# indices
begin
    iX = Vector{Int64}()
    iXglobal = Vector{Int64}()
    iY = Vector{Int64}()
    iYglobal = Vector{Int64}()
    iZ = Vector{Int64}()
    iZglobal = Vector{Int64}()

    iA = Vector{Int64}()
    iAglobal = Vector{Int64}()

    for (i, var) in enumerate(vars)
        if typeof(var) == SpatialVariable
            if var.axis == :X
                push!(iX, var.i)
                push!(iXglobal, i)
            elseif var.axis == :Y
                push!(iY, var.i)
                push!(iYglobal, i)
            else
                push!(iZ, var.i)
                push!(iZglobal, i)
            end
        else
            push!(iA, var.i)
            push!(iAglobal, i)
        end
    end
end

struct Indexer
    iX::Vector{Int64}
    iXg::Vector{Int64}
    iY::Vector{Int64}
    iYg::Vector{Int64}
    iZ::Vector{Int64}
    iZg::Vector{Int64}
    iA::Vector{Int64}
    iAg::Vector{Int64}
end

indexer = Indexer(iX, iXglobal, iY, iYglobal, iZ, iZglobal, iA, iAglobal)

function obj(variables::Vector{Float64}, p::TestParams, ind::Indexer)
    # update spatial variables
    Xnew = addvalues(p.X, ind.iX, variables[ind.iXg])
    Ynew = addvalues(p.Y, ind.iY, variables[ind.iYg])
    Znew = addvalues(p.Z, ind.iZ, variables[ind.iZg])

    # update internal variables
    Anew = replacevalues(p.A, ind.iA, variables[ind.iAg])

    # new stiffness matrices
    ks = [AsapToolkit.kglobal(Xnew, Ynew, Znew, e, a, id) for (e, a, id) in zip(p.E, Anew, p.nodeids)]

    K = assembleglobalK(ks, p)

    U = solveU(K, p)

    U' * p.P[p.freeids]
end

@time obj(initialvalues, p, indexer);

@time Zygote.gradient(var -> obj(var, p, indexer), initialvalues)[1]

objclosure(vals::Vector{Float64}, p::TestParams) = obj(vals, p, indexer);
@time objclosure(initialvalues, p)

func = Optimization.OptimizationFunction(objclosure, Optimization.AutoZygote())
prob = Optimization.OptimizationProblem(func, initialvalues, p;
    lb = lowerbounds,
    ub = upperbounds)


trace = Vector{Float64}()
hist = Vector{Vector{Float64}}()

function cb(vals::Vector{Float64}, loss::Float64)
    push!(trace, loss)
    push!(hist, deepcopy(vals))
    false
end

begin
    empty!(trace)
    empty!(hist)

    @time sol = Optimization.solve(prob,
        OptimizationNLopt.NLopt.LD_LBFGS(),
        reltol = 1e-4,
        callback = cb)
end

Xfinal = addvalues(p.X, indexer.iX, sol.u[indexer.iXg])
Yfinal = addvalues(p.Y, indexer.iY, sol.u[indexer.iYg])
Zfinal = addvalues(p.Z, indexer.iZ, sol.u[indexer.iZg])

Afinal = replacevalues(p.A, indexer.iA, sol.u[indexer.iAg])
Anormalized = Afinal ./ maximum(abs.(Afinal))
#plotting
begin
    p0 = Point3.(getproperty.(model.nodes, :position))
    e0 = vcat([p0[id] for id in p.nodeids]...)

    p1 = Point3.(Xfinal, Yfinal, Zfinal)
    e1 = vcat([p1[id] for id in p.nodeids]...)
end

begin
    fig = Figure()
    ax0 = Axis3(fig[1,1],
        aspect = :data)

    # hidedecorations!(ax0); hidespines!(ax0)

    linesegments!(e0, color = :white)

    ax1 = Axis3(fig[1,2],
        aspect = :data)

    hidedecorations!(ax1); hidespines!(ax1)

    linesegments!(e1, 
    color = (blue, 0.75),
    linewidth = Anormalized .* 10,
    )

    linkaxes!(ax0, ax1)

    fig
end