module  TestPCA

using PopGen
using Test

cats = @nancycats;

@testset "PCA.jl" begin
    @testset "pca" begin
        @test size(pca(cats).proj) == (237,93)
        @test size(pca(cats, maxpc = 25).proj) == (237,25)
        @test size(pca(cats, maxpc = 25).proj) == (237,25)
        @test size(pca(cats, maxpc = 25, center = true).proj) == (237,25)
        @test size(pca(cats, maxpc = 25, scale = false).proj) == (237,25)
        @test size(pca(cats, pratio = 0.7).proj) == (237,17)
    end
end

end #module