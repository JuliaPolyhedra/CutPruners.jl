################################################################################
# CutPruners
# A package to manage polyhedral convex functions
################################################################################

__precompile__()

module CutPruners

using MathProgBase
using DocStringExtensions

# Redudancy checking
include("redund.jl")

include("abstract.jl")
include("avg.jl")
include("decay.jl")
include("dematos.jl")
include("exact.jl")

end # module
