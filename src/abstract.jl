################################################################################
# Implement abstract type of CutPruner
################################################################################
export AbstractCutPruningAlgo
export CutPruner, AbstractCutPruner, addcuts!, ncuts

abstract AbstractCutPruningAlgo

"""
A cut pruner maintains a matrix `A` and a vector `b` such that
represents `size(A, 1)` (` == length(b)`) cuts.
Let `a_i` be `A[i,:]` and `β_i` be `b[i]`, the meaning of the cut depends on the sense.
Cuts (A, b) defines the half-space satisfying:
Ax >= b if islb
Ax <= b otherwise
If `sense` is
* `:Min`, then the cut pruner represents the concave polyhedral function `min ⟨a_i, x⟩ + β_i`;
* `:Max`, then the cut pruner represents the convex polyhedral function `max ⟨a_i, x⟩ + β_i`;
* `:≤`, then the cut pruner represents the polyhedra defined by the intersection of the half-space `⟨a_i, x⟩ ≤ β_i`;
* `:≥`, then the cut pruner represents the polyhedra defined by the intersection of the half-space `⟨a_i, x⟩ ≥ β_i`.

Internally, instead of `sense`, the booleans `isfun` and `islb` are stored.
The mapping between `sense` and these two booleans is given by the following table

| `sense` | `isfun` | `islb` |
| ------- | ------- | ------ |
| Min     | true    | false  |
| Max     | true    | true   |
| ≤       | false   | false  |
| ≥       | false   | true   |

"""
abstract AbstractCutPruner{N, T}

function gettype(sense::Symbol)
    if sense == :Min
        true, false
    elseif sense == :Max
        true, true
    elseif sense == :≤
        false, false
    elseif sense == :≥
        false, true
    else
        throw(ArgumentError("Invalid value `$sense' for sense. It should be :Min, :Max, :≤ or :≥."))
    end
end
function getsense(isfun::Bool, islb::Bool)
    if isfun
        islb ? :Max : :Min
    else
        islb ? :≥ : :≤
    end
end
getsense(pruner::AbstractCutPruner) = getsense(isfun(pruner), islb(pruner))
isfun(pruner::AbstractCutPruner) = pruner.isfun
islb(pruner::AbstractCutPruner) = pruner.islb

immutable CutPruner{N, T} end

"""Return whether the CutPruner `man` has any cut."""
function Base.isempty(man::AbstractCutPruner)
    return isempty(man.b)
end

"""Return number of cuts in CutPruner `man`."""
function ncuts(man::AbstractCutPruner)
    return length(man.b)
end

# COMPARISON
"""Get current `trust` of CutPruner `man`."""
gettrust(man::AbstractCutPruner) = man.trust

function _indmin(a::Vector, tiebreaker::Vector)
    imin = 1
    for i in 2:length(a)
        if a[i] < a[imin] || (a[i] == a[imin] && tiebreaker[i] < tiebreaker[imin])
            imin = i
        end
    end
    imin
end

"""Remove `num` cuts in CutPruner `man`."""
function choosecutstoremove(man::AbstractCutPruner, num::Int)
    # MergeSort is stable so in case of equality, the oldest cut loose
    # However PartialQuickSort is a lot faster

    trust = gettrust(man)
    if num == 1
        [_indmin(trust, man.ids)]
    else
        # /!\ PartialQuickSort is unstable, here it does not matter
        function _lt(i, j)
            # If cuts have same trust, remove oldest cut
            if trust[i] == trust[j]
                man.ids[i] < man.ids[j]
            # Else, remove cuts with lowest trust
            else
                trust[i] < trust[j]
            end
        end
        # Return index of `num` cuts with lowest trusts
        sort(1:length(trust), alg=PartialQuickSort(num), lt=_lt)[1:num]
    end
end

"""Test if cut `i` is better than `newcuttrust`."""
isbetter(man::AbstractCutPruner, i::Int, mycut::Bool) = gettrust(man)[i] > initialtrust(man, mycut)

# CHANGE

# Add cuts Ax >= b
# If mycut then the cut has been added because of one of my trials
function addcuts!{N, T}(man::AbstractCutPruner{N, T},
                     A::AbstractMatrix{T},
                     b::AbstractVector{T},
                     mycut::AbstractVector{Bool})
    # get current number of cuts:
    ncur = ncuts(man)
    nincumbents = size(A, 1)

    # check redundancy
    redundants = checkredundancy(man.A, man.b, A, b, man.isfun, man.islb, man.TOL_EPS)
    tokeep = setdiff(collect(1:nincumbents), redundants)

    # if all cuts are redundants, then do nothing:
    if length(tokeep) == 0
        return zeros(Int, nincumbents)
    end
    A = A[tokeep, :]
    b = b[tokeep]

    # get number of new cuts in A:
    nnew = size(A, 1)
    @assert length(mycut) == length(b) == nnew
    @assert nnew > 0

    # If not enough room, need to remove some cuts
    if man.maxncuts != -1 && ncur + nnew > man.maxncuts
        # get indexes of cuts with lowest trusts:
        R = choosecutstoremove(man, ncur + nnew - man.maxncuts)

        # Check if some new cuts should be ignored
        take = man.maxncuts - ncur
        # get number of cuts to remove
        nreplaced = length(R)
        nmycut = sum(mycut)
        # Start:
        # |      | J      j|
        # | take |         |
        # |  mycut  |      |
        # End:
        # |      | J   j   |
        # | take       |   |
        # |  mycut  |      |
        while take + length(R) - nreplaced < nnew
            # I first try to see if the nmycut cuts generated by me can be taken
            # because they have better trust than the nnew-nmycut others
            if isbetter(man, R[nreplaced], take < nmycut)
                nreplaced -= 1
            else
                take += 1
            end
        end
        # Cuts that will be replaced
        R = @view R[1:nreplaced]
        # Nowe we split A, b into
        # * A, b   : cuts to be pushed
        # * Ar, br : cuts to replace old cuts
        # * _, _   : cuts removed
        if nreplaced == nnew
            status = R
            Ar = A
            br = b
            mycutr = mycut
            A = similar(A, 0, size(A, 2))
            b = similar(b, 0)
            mycut = similar(mycut, 0)
        else
            status = zeros(Int, nnew)
            if take < nnew
                # Remove ignored cuts
                takemycut = min(take, nmycut)
                takenotmycut = take - takemycut
                takeit = zeros(Bool, nnew)
                for i in 1:nnew
                    ok = false
                    if mycut[i]
                        if takemycut > 0
                            takemycut -= 1
                            ok = true
                        end
                    else
                        if takenotmycut > 0
                            takenotmycut -= 1
                            ok = true
                        end
                    end
                    if ok
                        takeit[i] = true
                    end
                end

                takeit = find(takeit)
            else
                takeit = collect(1:nnew)
            end
            @assert take == length(takeit)
            replaced = takeit[1:nreplaced]
            pushed = takeit[(nreplaced+1):end]

            status[replaced] = R
            Ar = @view A[replaced,:]
            br = @view b[replaced]
            mycutr = @view mycut[replaced]
            status[pushed] = ncur + (1:length(pushed))
            A = @view A[pushed,:]
            b = @view b[pushed]
            mycut = @view mycut[pushed]
        end
        if nreplaced > 0
            man.A[R, :] = Ar
            man.b[R] = br
            replacecuts!(man, R, mycutr)
        end
    end
    if !isempty(b)
        # Just append cuts
        status = ncur + (1:nnew)
        man.A = [man.A; A]
        man.b = [man.b; b]
        pushcuts!(man, mycut)
    end

    status
end

"""Keep only cuts whose indexes are in Vector `K`."""
function keeponly!(man::AbstractCutPruner, K::AbstractVector{Int})
    man.A = man.A[K, :]
    man.b = man.b[K]
    man.trust = man.trust[K]
end

"""Get a Vector of Float64 specifying the initial trusts of `mycut`."""
function initialtrusts(man::AbstractCutPruner, mycut::AbstractVector{Bool})
    Float64[initialtrust(man, mc) for mc in mycut]
end

"""Reset trust of cuts with indexes in `js`."""
function replacecuts!(man::AbstractCutPruner, js::AbstractVector{Int}, mycut::AbstractVector{Bool})
    man.trust[js] = initialtrusts(man, mycut)
    man.ids[js] = newids(man, length(js))
end

function pushcuts!(man::AbstractCutPruner, mycut::AbstractVector{Bool})
    append!(man.trust, initialtrusts(man, mycut))
    append!(man.ids, newids(man, length(mycut)))
end

function newids(man::AbstractCutPruner, n::Int)
    (man.id+1):(man.id += n)
end

function normalizedcut{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, k::Int, isfun::Bool, tol::Float64)
    a = @view A[k, :]
    β = b[k]
    if isfun || abs(β) < tol
        a, β
    else
        a / β, one(T)
    end
end

"""
Check redundant cuts. Return index of redundant cuts in `Anew`.

$(SIGNATURES)

"""
function checkredundancy{T}(A::AbstractMatrix{T}, b::AbstractVector{T},
                            Anew::AbstractMatrix{T}, bnew::AbstractVector{T},
                            isfun::Bool, islb::Bool, tol::Float64)
    # index of redundants cuts
    redundants = Int[]
    # number of new lines
    nnew = size(Anew, 1)

    for kk in 1:nnew
        a, β = normalizedcut(Anew, bnew, kk, isfun, tol)
        chk, indk = isinside(A, b, a, isfun, tol)
        if chk
            ared, βred = normalizedcut(A, b, indk, isfun, tol)
            if islb ? β <= βred+tol : β+tol >= βred
                push!(redundants, kk)
            end
        end
    end

    redundants
end


"""Check if `λ` is a line of matrix `A`. `λ` might not have the same `eltype` as `A` and `b` as it might have been scaled by `normalizecut`."""
function isinside{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, λ::AbstractVector, isfun::Bool, tol::Float64)
    nlines = size(A, 1)

    check = false
    k = 0
    while ~check && k < nlines
        k += 1
        a, β = normalizedcut(A, b, k, isfun, tol)
        check = norm(a - λ, Inf) < tol
    end
    check, k
end
