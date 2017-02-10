################################################################################
# CutPruner
# A package to manage polyhedral convex functions
################################################################################

__precompile__()

module CutPruner

using DocStringExtensions

include("abstract.jl")
include("avg.jl")
include("decay.jl")

end # module
