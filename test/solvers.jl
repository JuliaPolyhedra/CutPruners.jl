# Inspired from JuMP/test/solvers

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

# Load available solvers
grb = try_import(:Gurobi)
cpx = try_import(:CPLEX)
xpr = try_import(:Xpress)
mos = false && try_import(:MosekTools)
cbc = try_import(:Cbc)
if cbc; import Clp; end

import JuMP
import GLPK
glp = true
ipt = try_import(:Ipopt)
eco = try_import(:ECOS)

# Create LP solver list
lp_solvers = Module[]
grb && push!(lp_solvers, Gurobi)
cpx && push!(lp_solvers, CPLEX)
xpr && push!(lp_solvers, Xpress)
mos && push!(lp_solvers, Mosek)
cbc && push!(lp_solvers, Clp)
glp && push!(lp_solvers, GLPK)
ipt && push!(lp_solvers, Ipopt)
eco && push!(lp_solvers, ECOS)
