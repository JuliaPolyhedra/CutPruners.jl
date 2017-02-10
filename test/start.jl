@testset "Initialization" begin
    for pruner in [AvgCutPruner{Int}(1), DecayCutPruner{Int}(1)]
        @test !isstarted(pruner)
        start!(pruner, 0)
        @test ncuts(pruner) == 0
        @test isstarted(pruner)
    end
end
