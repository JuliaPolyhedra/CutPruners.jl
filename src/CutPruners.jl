################################################################################
# CutPruners
# A package to manage polyhedral convex functions
################################################################################

__precompile__()

module CutPruners

using JuMP, DocStringExtensions

include("abstract.jl")
include("avg.jl")
include("decay.jl")
include("level1.jl")
include("exact.jl")

end # module
