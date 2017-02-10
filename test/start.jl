@testset "Initialization" begin
    for algo in [AvgCutPruningAlgo(1), DecayCutPruningAlgo(1)]
        pruner = CutPruner{2, Int}(algo)
        @test ncuts(pruner) == 0
    end
end
