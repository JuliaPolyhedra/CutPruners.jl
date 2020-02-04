################################################################################
# CutPruners
# A package to manage polyhedral convex functions
################################################################################


module CutPruners

using Compat, Compat.LinearAlgebra, Compat.SparseArrays
using JuMP

# Redudancy checking
include("redund.jl")

include("abstract.jl")
include("avg.jl")
include("decay.jl")
include("dematos.jl")
include("exact.jl")

end # module
