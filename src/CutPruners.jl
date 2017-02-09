################################################################################
# CutPruners
# A package to manage polyhedral convex functions
################################################################################

__precompile__()

module CutPruners

using DocStringExtensions

include("abstract.jl")
include("avg.jl")
include("decay.jl")
include("level1.jl")

end # module
