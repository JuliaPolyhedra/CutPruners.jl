using CutPruners
using Base.Test

@testset "Level-1 cut pruning" begin
    algo = LevelOnePruningAlgo(3)
    pruner = CutPruner{2, Int}(algo)
    addcuts!(pruner, [1 0], [0], [true])
    addcuts!(pruner, [2 0], [0], [true])
    addcuts!(pruner, [3 0], [0], [true])
    # FIXME: need to update stats before reaching the ncuts limit ...
    CutPruners.updatestats!(pruner, [1 0])
    addcuts!(pruner, [4 0], [0], [true])
    @test sort(pruner.cuts_DE[:,1]) == [2, 3, 4]
    CutPruners.updatestats!(pruner, [1 0])
    addcuts!(pruner, [5 0], [0], [true])
    @test sort(pruner.cuts_DE[:,1]) == [3, 4, 5]

    @test ncuts(pruner) == 3

    @test ~CutPruners.checkconsistency(pruner)
    CutPruners.updatestats!(pruner, [1 0])
    @test CutPruners.checkconsistency(pruner)
end
