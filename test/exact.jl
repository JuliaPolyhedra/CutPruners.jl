@testset "Exact cuts pruning" begin
    algo = ExactPruningAlgo(3, MathProgBase.DefaultLPSolver(),
                            [-10, -10], [10, 10])
    pruner = CutPruner{2, Int}(algo)
    addcuts!(pruner, [1 0], [0], [true])
    addcuts!(pruner, [2 0], [0], [true])
    addcuts!(pruner, [3 0], [0], [true])
    addcuts!(pruner, [4 0], [0], [true])
    @test sort(pruner.cuts_DE[:,1]) == [2, 3, 4]
    addcuts!(pruner, [5 0], [0], [true])
    @test sort(pruner.cuts_DE[:,1]) == [3, 4, 5]
    # However, order should only be used to break ties
    addcuts!(pruner, [6 0], [0], [false])
    @test sort(pruner.cuts_DE[:,1]) == [4, 5, 6]
    addcuts!(pruner, [7 0; 8 0], [0, 0], [true, true])
    @test sort(pruner.cuts_DE[:,1]) == [6, 7, 8]
end
