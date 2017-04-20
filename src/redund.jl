function normalizedcut{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, k::Int, isfun::Bool, tol::Float64)
    a = @view A[k, :]
    β = b[k]
    na = norm(a, 2)
    if isfun || na < tol
        a, β
    else
        a / na, β / na
    end
end

"""
Check redundant cuts. Return index of redundant cuts in `Anew`.

$(SIGNATURES)

"""
function checkredundancy{T}(A::AbstractMatrix{T}, b::AbstractVector{T},
                            Anew::AbstractMatrix{T}, bnew::AbstractVector{T},
                            isfun::Bool, islb::Bool, tol::Float64, ident::Bool=false)
    # index of redundants cuts
    redundants = Int[]
    # number of new lines
    nnew = size(Anew, 1)

    for kk in 1:nnew
        a, β = normalizedcut(Anew, bnew, kk, isfun, tol)
        chk, indk = isinside(A, b, a, isfun, tol)
        if chk && (~ident || indk!=kk)
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
        check = (norm(a - λ, Inf) < tol)
    end
    check, k
end


checkredundancy(A, b, isfun, islb, tol)=checkredundancy(A, b, A, b, isfun, islb, tol, true)


"""Remove redundants cuts in polyhedra (A, b)."""
function clean!(A, b, isfun, islb, tol)
    nincumbents = size(A, 1)
    redundants = checkredundancy(A, b, isfun, islb, tol)
    if !isempty(redundants) && length(redundants) < nincumbents
        tokeep = setdiff(collect(1:nincumbents), redundants)

        A = A[tokeep, :]
        b = b[tokeep]
    end
end

