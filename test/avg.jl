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
