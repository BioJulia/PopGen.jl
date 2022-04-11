module  TestKMeans

using PopGen
using Test

cats = @nancycats;

@testset "KMeans.jl" begin
    @testset "kmeans" begin
    @test kmeans(cats, 2:5) isa PopGen.KMeansResults
    @test kmeans(cats, [2,5,6]) isa PopGen.KMeansResults
    @test kmeans(cats, 2:5, 30) isa PopGen.KMeansResults
    @test kmeans(cats, [3,6,7], 30) isa PopGen.KMeansResults
    @test kmeans(cats, krange = 2:5, iterations = 30) isa PopGen.KMeansResults
    @test kmeans(cats, krange = [2,3,4], iterations = 30) isa PopGen.KMeansResults
    tmp = kmeans(cats, krange = 2:5, iterations = 30)
    @test size(tmp.assignments) == (237,5)
    end
end

end #module