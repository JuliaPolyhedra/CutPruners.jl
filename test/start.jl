@testset "Initialization" begin
    for algo in [AvgCutPruningAlgo(1), DecayCutPruningAlgo(1), LevelOnePruningAlgo(3)]
        @test_throws ArgumentError CutPruner{2, Int}(algo, :Mix)
        for sense in [:Min, :Max, :≤, :≥]
            pruner = CutPruner{2, Int}(algo, sense)
            @test ncuts(pruner) == 0
            @test CutPruners.isfun(pruner) == (sense in [:Min, :Max])
            @test CutPruners.islb(pruner) == (sense in [:Max, :≥])
            @test CutPruners.getsense(pruner) == sense
        end
    end
end

@testset "Check redundancy" begin
    # check redundant
    TOL = 1e-6
    A = rand(10, 2)
    b = rand(10)
    ck = CutPruners.isinside(A, b, A[1,:], true, TOL)
    @test ck[1]

    Anew = A[9:10, :]
    bnew = b[9:10]
    red = CutPruners.checkredundancy(A, b, Anew, bnew, true, true, TOL)
    @test red == [1, 2]
    A = rand(10, 2)
    ck = CutPruners.isinside(A, b, A[1,:], true, TOL)
    @test ck[1]

    # check non redundant
    b = rand(10)
    Anew = rand(2, 2)
    ck = CutPruners.isinside(A, b, Anew[1, :], true, TOL)
    @test ~ck[1]
    bnew = rand(2)
    red = CutPruners.checkredundancy(A, b, Anew, bnew, true, true, TOL)
    @test red == []
end
