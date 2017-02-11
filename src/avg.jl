export AvgCutPruningAlgo, AvgCutPruner

"""
$(TYPEDEF)

Removes the cuts with lower trust where the trust is: nused / nwith + bonus
where the cut has been used `nused` times amoung `nwith` optimization done with it.
We say that the cut was used if its dual value is nonzero.
It has a bonus equal to `mycutbonus` if the cut was generated using a trial given by the problem using this cut.
If `nwidth` is zero, `nused/nwith` is replaced by `newcuttrust`.
"""
type AvgCutPruningAlgo <: AbstractCutPruningAlgo
    # maximum number of cuts
    maxncuts::Int
    newcuttrust::Float64
    mycutbonus::Float64
    function AvgCutPruningAlgo(maxncuts::Int, newcuttrust=3/4, mycutbonus=1/4)
        new(maxncuts, newcuttrust, mycutbonus)
    end
end

type AvgCutPruner{N, T} <: AbstractCutPruner{N, T}
    # used to generate cuts
    # Cuts (A, b) defines the half-space satisfying: Ax >= b
    # A
    cuts_DE::AbstractMatrix{T}
    # b
    cuts_de::AbstractVector{T}

    # number of optimality cuts
    nσ::Int
    # number of feasibility cuts
    nρ::Int
    # index of optimality cuts
    σs::Vector{Int}
    # index of feasibility cuts
    ρs::Vector{Int}

    # maximum number of cuts
    maxncuts::Int

    # number of optimization performed
    nwith::Vector{Int}
    # number of times where the cuts have been used
    nused::Vector{Int}
    mycut::Vector{Bool}
    trust::Nullable{Vector{Float64}}
    ids::Vector{Int} # small id means old
    id::Int # current id

    newcuttrust::Float64
    mycutbonus::Float64

    function AvgCutPruner(maxncuts::Int, newcuttrust=3/4, mycutbonus=1/4)
        new(spzeros(T, 0, N), T[], 0, 0, Int[], Int[], maxncuts, Int[], Int[], Bool[], nothing, Int[], 0, newcuttrust, mycutbonus)
    end
end

(::Type{CutPruner{N, T}}){N, T}(algo::AvgCutPruningAlgo) = AvgCutPruner{N, T}(algo.maxncuts, algo.newcuttrust, algo.mycutbonus)

# COMPARISON
"""Update cuts relevantness after a solver's call returning dual vector `σ\rho`."""
function updatestats!(man::AvgCutPruner, σρ)
    if ncuts(man) > 0
        man.nwith += 1
        # TODO: dry 1e-6 in CutPruner?
        man.nused[σρ .> 1e-6] += 1
        man.trust = nothing # need to be recomputed
    end
end

function gettrustof(man::AvgCutPruner, nwith, nused, mycut)
    (nwith == 0 ? man.newcuttrust : nused / nwith) + (mycut ? man.mycutbonus : 0)
end
function initialtrust(man::AvgCutPruner, mycut)
    gettrustof(man, 0, 0, mycut)
end
function gettrust(man::AvgCutPruner)
    if isnull(man.trust)
        trust = man.nused ./ man.nwith
        trust[man.nwith .== 0] = man.newcuttrust
        trust[man.mycut] += man.mycutbonus
        man.trust = trust
    end
    get(man.trust)
end

# CHANGE

#FIXME: do not drop cuts in cuts_DE and cuts_de?
function keeponly!(man::AvgCutPruner, K::AbstractVector{Int})
    man.nwith = man.nwith[K]
    man.nused = man.nused[K]
    man.mycut = man.mycut[K]
    man.trust = gettrust(man)[K]
end

function replacecuts!(man::AvgCutPruner, js::AbstractVector{Int}, mycut::AbstractVector{Bool})
    man.nwith[js] = 0
    man.nused[js] = 0
    man.mycut[js] = mycut
    gettrust(man)[js] = initialtrusts(man, mycut)
    man.ids[js] = newids(man, length(js))
end

"""Push new cut in CutPruner `man`."""
function pushcuts!(man::AvgCutPruner, mycut::AbstractVector{Bool})
    n = length(mycut)
    append!(man.nwith, zeros(n))
    append!(man.nused, zeros(n))
    append!(man.mycut, mycut)
    if !isnull(man.trust)
        append!(get(man.trust), initialtrusts(man, mycut))
    end
    append!(man.ids, newids(man, n))
end
