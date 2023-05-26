@reexport using ChainRulesCore
@reexport import Optimization, OptimizationNLopt
@reexport using IterativeSolvers

include("Truss/Translation.jl")

include("Truss/Types.jl")
export TrussOptProblem
export SpatialVariable, AreaVariable, CoupledVariable, AbstractVariable

include("Truss/Functions.jl")
export kglobal, L, Rtruss, assembleglobalK, solveU, Utruss, replacevalues, addvalues

include("Truss/Adjoints.jl")