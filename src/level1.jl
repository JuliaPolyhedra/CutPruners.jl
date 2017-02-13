export LevelOneCutPruner, LevelOnePruningAlgo


type LevelOnePruningAlgo <: AbstractCutPruningAlgo
    # maximum number of cuts
    maxncuts::Int
    function LevelOnePruningAlgo(maxncuts::Int)
        new(maxncuts)
    end
end


"""
$(TYPEDEF)

Removes the cuts with lower trust where the trust is: nused / nwith + bonus
where the cut has been used `nused` times amoung `nwith` optimization done with it.
We say that the cut was used if its dual value is nonzero.
It has a bonus equal to `mycutbonus` if the cut was generated using a trial given by the problem using this cut.
If `nwidth` is zero, `nused/nwith` is replaced by `newcuttrust`.
"""
type LevelOneCutPruner{N, T} <: AbstractCutPruner{N, T}
    # used to generate cuts
    cuts_DE::AbstractMatrix{T}
    cuts_de::AbstractVector{T}

    maxncuts::Int

    trust::Nullable{Vector{Float64}}
    ids::Vector{Int} # small id means old
    id::Int # current id

    #set of states where cut k is active
    territories::Vector{Vector{Tuple{Int64, T}}}
    nstates::Int
    states::Array{T, 2}

    function LevelOneCutPruner(maxncuts::Int)
        new(spzeros(T, 0, N), T[], maxncuts, nothing, Int[], 0, [], 0, zeros(T, 0, N))
    end
end

(::Type{CutPruner{N, T}}){N, T}(algo::LevelOnePruningAlgo) = LevelOneCutPruner{N, T}(algo.maxncuts)


"""Update territories with cuts previously computed during backward pass.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`
* `position::Array{T, 2}`
    New visited positions
"""
function updatestats!{T}(man::LevelOneCutPruner, position::Array{T, 2})
    nc = ncuts(man)
    # get number of new positions to analyse:
    nx = size(position, 1)
    # get number of cuts added since last update:
    nt = nc - length(man.territories)

    man.territories = vcat(man.territories, [Tuple{Int64, T}[] for _ in 1:nt])
    for i in 1:nt
        # update territory for new cut with index i
        updateterritory!(man, nc - nt + i)
    end

    for i in 1:nx
        addstate!(man, position[i, :])
    end

    # need to be recomputed
    man.trust = nothing
end


"""Add a new state to test and accordingly update territories of each cut.

$(SIGNATURES)

"""
function addstate!(man::LevelOneCutPruner, x::Vector)
    # Get cut which is active at point `x`:
    bcost, bcuts = optimalcut(man, x)

    # update number of states
    man.nstates += 1
    # Add `x` with index nstates  to the territory of cut with index `bcuts`:
    push!(man.territories[bcuts], (man.nstates, bcost))

    # Add `x` to the list of visited state:
    man.states = vcat(man.states, x')
end


"""Find active cut at point `xf`.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`:
    CutPruner
* `xf::Vector{Float64}`:

# Return
`bestcost::Float64`
    Value of supporting cut at point `xf`
`bestcut::Int64`
    Index of supporting cut at point `xf`
"""
function optimalcut{T}(man::LevelOneCutPruner,
                       xf::Vector{T})
    bestcost = -Inf::Float64
    bestcut = -1
    dimstates = length(xf)
    nc = ncuts(man)

    @inbounds for i in 1:nc
        cost = -man.cuts_de[i]
        for j in 1:dimstates
            cost += man.cuts_DE[i, j]*xf[j]
        end
        if cost > bestcost
            bestcost = cost
            bestcut = i
        end
    end
    return bestcost, bestcut
end


"""Update territories (i.e. the set of tested states where
    a given cut is active) considering new cut given by index `indcut`.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`:
* `indcut::Int64`:
    new cut index
"""
function updateterritory!(man::LevelOneCutPruner, indcut::Int64)
    for k in 1:ncuts(man)
        if k == indcut
            continue
        end
        todelete = []
        for (num, (ix, cost)) in enumerate(man.territories[k])
            x = man.states[ix, :]

            costnewcut = cutvalue(man, indcut, x)

            if costnewcut > cost
                push!(todelete, num)
                push!(man.territories[indcut], (ix, costnewcut))
            end
        end
        deleteat!(man.territories[k], todelete)
    end
end


"""
Get value of cut with index `indc` at point `x`.

$(SIGNATURES)

# Arguments
- `man::LevelOneCutPruner`
    Approximation of the value function as linear cuts
- `indc::Int64`
    Index of cut to consider
- `x::Array{Float64}`
    Coordinates of state

# Return
`cost::Float64`
    Value of cut `indc` at point `x`
"""
function cutvalue(man::LevelOneCutPruner, indc::Int, x::Vector{Float64})
    cost = -man.cuts_de[indc]
    for j in 1:length(x)
        cost += man.cuts_DE[indc, j]*x[j]
    end
    cost
end


function gettrustof(man::LevelOneCutPruner, nwith, nused, mycut)
    (nwith == 0 ? man.newcuttrust : nused / nwith) + (mycut ? man.mycutbonus : 0)
end

initialtrust(man::LevelOneCutPruner, mycut) = 0

function gettrust(man::LevelOneCutPruner)
    if isnull(man.trust) || ~checkconsistency(man)
        trust = Float64[length(terr) for terr in man.territories]
        man.trust = trust
    end
    get(man.trust)
end


function keeponly!(man::LevelOneCutPruner, K::AbstractVector{Int})
    man.trust = gettrust(man)[K]
end


function replacecuts!(man::LevelOneCutPruner, js::AbstractVector{Int}, mycut::AbstractVector{Bool})
    gettrust(man)[js] = initialtrusts(man, mycut)
    man.territories = man.territories[js]
    man.ids[js] = newids(man, length(js))
end


"""Push new cut in CutPruner `man`."""
function pushcuts!(man::LevelOneCutPruner, mycut::AbstractVector{Bool})
    n = length(mycut)
    if !isnull(man.trust)
        append!(get(man.trust), initialtrusts(man, mycut))
    end
    append!(man.ids, newids(man, n))
end


"""Check consistency of CutPruner `man`."""
function checkconsistency(man::LevelOneCutPruner)
    consistency = (ncuts(man) == length(man.territories))
    consistency
end
