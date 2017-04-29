################################################################################
# CutPruners
# A package to manage polyhedral convex functions
################################################################################

__precompile__()

module CutPruners

using Compat
using DocStringExtensions

# Redudancy checking
include("redund.jl")

include("abstract.jl")
include("avg.jl")
include("decay.jl")
include("dematos.jl")

end # module
