
export exactpruning!

# Exact pruning
"""Remove dominated cuts in CutPruner `man`.

We use a LP solver to determine whether a cut is dominated or not.
"""
function exactpruning!(man::AbstractCutPruner, solver::MathProgBase.AbstractMathProgSolver)
    K = getdominated(man.A, man.b, man.islb, man.isfun, solver)
    removecuts!(man, K)
end

"""Return dominated cuts."""
function getdominated(A, b, islb, isfun, solver)
    red = Int[]
    if size(A, 1) == 1
        return red
    end
    for i in 1:size(A, 1)
        if isdominated(A, b, islb, isfun, i, solver)
            push!(red, i)
        end
    end
    red
end

"""State whether a cut is dominated with a tolerance epsilon."""
function isdominated(A, b, islb, isfun, k, solver, epsilon=1e-8)
    # we use MathProgBase to solve the test
    # The LP is:
    # min - y
    #     x ∈ R^n, y ∈ R
    #     y <= (A[k, :] - A[j, :]) + b[k] - b[j]     ∀ i != k

    # we formulate an equivalent problem as
    # min c'*z
    # s.t H*z <= h

    # get problem's dimension
    ncuts, nx = size(A)
    # allocate arrays
    h = zeros(ncuts - 1)
    H = zeros(ncuts - 1, nx + 1)
    c = zeros(nx + 1)
    c[1] = -1

    λk = @view A[k, :]

    ic = 0
    @inbounds for ix in 1:ncuts
        if ix != k
            ic += 1
            h[ic] = b[k] - b[ix]
            H[ic, 1] = 1

            @inbounds for jx in 1:nx
                H[ic, jx+1] = A[ic, jx] - λk[jx]
            end
        end
    end

    # solve the LP with MathProgBase
    res = linprog(c, H, -Inf, h,  -Inf, Inf, solver).objval
    return (-res < epsilon)
end
