@testset "Initialization" begin
    for algo in [AvgCutPruningAlgo(1), DecayCutPruningAlgo(1)]
        pruner = CutPruner{2, Int}(algo)
        @test ncuts(pruner) == 0
    end
end

@testset "Check redundancy" begin
    # check redundant
    TOL = 1e-6
    A = rand(10, 2)
    ck = CutPruners.isinside(A, A[1,:], TOL)
    @test ck[1]

    b = rand(10)
    Anew = A[9:10, :]
    bnew = b[9:10]
    red = CutPruners.checkredundancy(A, b, Anew, bnew, TOL)
    @test red == [1, 2]
    A = rand(10, 2)
    ck = CutPruners.isinside(A, A[1,:], TOL)
    @test ck[1]

    # check non redundant
    b = rand(10)
    Anew = rand(2, 2)
    ck = CutPruners.isinside(A, Anew[1, :], TOL)
    @test ~ck[1]
    bnew = rand(2)
    red = CutPruners.checkredundancy(A, b, Anew, bnew, TOL)
    @test red == []
end
