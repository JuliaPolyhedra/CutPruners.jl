export ExactCutPruner, ExactPruningAlgo


type ExactPruningAlgo <: AbstractCutPruningAlgo
    # maximum number of cuts
    maxncuts::Int
    # LP solver to find dominant cut
    solver
    # lower bound of primal space
    lbounds::Vector
    # upper bound of primal space
    ubounds::Vector
    # solver's tolerance
    tol::Float64

    function ExactPruningAlgo(maxncuts::Int, solver,
                              lbounds, ubounds;
                              tolerance::Float64=1e-5)
        new(maxncuts, solver, lbounds, ubounds, tolerance)
    end
end


"""
$(TYPEDEF)

"""
type ExactCutPruner{N, T} <: AbstractCutPruner{N, T}
    # used to generate cuts
    cuts_DE::AbstractMatrix{T}
    cuts_de::AbstractVector{T}

    maxncuts::Int

    trust::Nullable{Vector{Bool}}
    ids::Vector{Int} # small id means old
    id::Int # current id

    solver
    lbounds::Vector{T}
    ubounds::Vector{T}
    TOL::Float64

    function ExactCutPruner(maxncuts::Int, solver, lbounds, ubounds, tol)
        new(spzeros(T, 0, N), T[], maxncuts, nothing, Int[], 0, solver, lbounds, ubounds, tol)
    end
end

(::Type{CutPruner{N, T}}){N, T}(algo::ExactPruningAlgo) = ExactCutPruner{N, T}(algo.maxncuts, algo.solver, algo.lbounds, algo.ubounds, algo.tol)


"""
Test whether the cut number k is dominated in CutPruner `man`.

$(SIGNATURES)

# Arguments
* `man::ExactCutPruner`:
* `k::Int`:
    Position of cut to test in CutPruner object

# Return
* `Bool`: true if the cut is dominant, false otherwise
"""
function isdominant(man::ExactCutPruner, k::Int)
    nx = length(man.lbounds)
    m = Model(solver=man.solver)
    @variable(m, alpha)
    @variable(m, man.lbounds[i] <= x[i=1:nx] <= man.ubounds[i])

    for i in 1:ncuts(man)
        if i!=k
            lambda = @view man.cuts_DE[i, :]
            @constraint(m, -man.cuts_de[i] + dot(lambda, x) <= alpha)
        end
    end

    λ_k = @view man.cuts_DE[k, :]
    @objective(m, Min, alpha - dot(λ_k, x) + man.cuts_de[k])
    solve(m)
    sol = getobjectivevalue(m)
    return (sol < man.TOL)
end


initialtrust(man::ExactCutPruner, mycut) = true

function gettrust(man::ExactCutPruner)
    if isnull(man.trust)
        trust = Bool[isdominant(man, i) for i in 1:ncuts(man)]
        man.trust = trust
    end
    get(man.trust)
end

# recompute trust only in positition specified in `pos`
# Use this function to not re-test all cuts but only a smaller subset
function gettrust(man::ExactCutPruner, pos::Vector{Int})
    if isnull(man.trust)
        trust = ones(Bool, ncuts(man))
        for pp in pos
            trust[pp] = isdominant(man, pp)
        end
        man.trust = trust
    end
    get(man.trust)
end


function keeponly!(man::ExactCutPruner, K::AbstractVector{Int})
    man.trust = gettrust(man)[K]
end


function replacecuts!(man::ExactCutPruner, js::AbstractVector{Int}, mycut::AbstractVector{Bool})
    gettrust(man)[js] = initialtrusts(man, mycut)
    man.ids[js] = newids(man, length(js))
end


"""Push new cut in CutPruner `man`."""
function pushcuts!(man::ExactCutPruner, mycut::AbstractVector{Bool})
    n = length(mycut)
    if !isnull(man.trust)
        append!(get(man.trust), initialtrusts(man, mycut))
    end
    append!(man.ids, newids(man, n))
end


"""Check consistency of CutPruner `man`."""
function checkconsistency(man::ExactCutPruner)
    # TODO: broken function
    consistency = (ncuts(man) == length(man.trust))
    consistency
end
