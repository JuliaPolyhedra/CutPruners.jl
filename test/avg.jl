@testset "Avg cut pruning" begin
    algo = AvgCutPruningAlgo(2)
    pruner = CutPruner{2, Int}(algo, :≤)
    addcuts!(pruner, [1 0], [1], [true])
    addcuts!(pruner, [0 1], [1], [true])
    CutPruners.updatestats!(pruner, [1, 0])
    addcuts!(pruner, [1 1; -1 -1; 0 1], [1, 1, 2], [true, false, true])
    @test pruner.A == [1 0; 1 1]
    @test pruner.b == [1, 1]
    @test pruner.ids == [1, 3]
end

@testset "Manual Cut Pruning" begin
    algo = AvgCutPruningAlgo(10)
    pruner = CutPruner{2, Int}(algo, :≤)
    addcuts!(pruner, [1 0], [1], [true])
    addcuts!(pruner, [0 1], [1], [true])
    CutPruners.updatestats!(pruner, [1 0])

    @test ncuts(pruner) == 2
    # change tolerance of pruning
    pruner.TOL_PRUNING = .5
    # prune cuts
    prunecuts!(pruner)
    @test ncuts(pruner) == 1
    # The cut [2 0] is taken even if [1 0] has one territory
    A, b = fetchcuts(pruner)
    @test A == [1 0]
    @test b == [1]
end
