@testset "Initialization" begin
    for man in [AvgCutPruner{Int}(1), DecayCutPruner{Int}(1)]
        @test !isstarted(man)
        start!(man, 0)
        @test isstarted(man)
    end
end
