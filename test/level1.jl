using CutPruners
using Base.Test

@testset "Level-1 cut pruning" begin
    algo = LevelOnePruningAlgo(3)
    pruner = CutPruner{2, Int}(algo, :Max)
    addcuts!(pruner, [1 0], [0], [true])
    addcuts!(pruner, [2 0], [0], [true])
    addcuts!(pruner, [3 0], [0], [true])
    # FIXME: need to update stats before reaching the ncuts limit ...
    CutPruners.updatestats!(pruner, [1 0])
    @test pruner.trust == [0, 0, 1.0]
    addcuts!(pruner, [4 0], [0], [true])
    @test pruner.A[:,1] == [4, 2, 3]
    @test pruner.trust == [1.0, 0, 0]
    CutPruners.updatestats!(pruner, [1 0])
    @test pruner.trust == [2.0, 0, 0]
    addcuts!(pruner, [5 0], [0], [true])
    @test pruner.A[:,1] == [4, 5, 3]
    @test pruner.trust == [0, 2.0, 0]

    @test ncuts(pruner) == 3

    CutPruners.updatestats!(pruner, [1 0])
    @test pruner.trust == [0, 3.0, 0]
end

@testset "Priority to new cuts" begin
    algo = LevelOnePruningAlgo(1)
    pruner = CutPruner{2, Int}(algo, :Min)
    addcuts!(pruner, [1 0], [0], [true])
    CutPruners.updatestats!(pruner, [1 0])
    addcuts!(pruner, [2 0], [0], [true])

    @test ncuts(pruner) == 1
    # The cut [2 0] is taken even if [1 0] has one territory
    @test pruner.A == [2 0]
end
