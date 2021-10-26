module  TestHeterozygosity

using PopGen
using Test

cats = @nancycats;
testarray = cats.genodata.genotype[1:10]

@testset "Heterozygosity.jl" begin
    @testset "counthet" begin
    @test PopGen.counthet(testarray, 135) == 5
    @test PopGen.counthet(testarray, 143) == 4
    @test PopGen.counthet(testarray, [135,143]) == [5,4]
    end

    @testset "counthom" begin
        @test PopGen.counthom(testarray, 135) == 2
        @test PopGen.counthom(testarray, 143) == 0
        @test PopGen.counthom(testarray, [135,143]) == [2,0]
    end
        
    @testset "het obs/exp" begin
        @test PopGen._hetero_obs(testarray) == 0.6
        @test PopGen._hetero_exp(testarray) == 0.6015625
    end

    @testset "Nei het" begin
        @test ismissing(PopGen._genediversitynei87(missing,0.5,11))
        @test ismissing(PopGen._genediversitynei87(0.5,missing,11))
        @test ismissing(PopGen._genediversitynei87(missing,missing,11))
        @test PopGen._genediversitynei87(0.5,0.5,11) == 0.525
        @test PopGen._genediversitynei87(0.5,0.5,11, corr = false) == 0.4772727272727273
    end

    @testset "heterozygosity" begin
        tmp = heterozygosity(cats)
        @test tmp isa DataFrames.DataFrame
        @test size(tmp) == (9,4)
        tmp = heterozygosity(cats, by = "sample")
        @test tmp isa DataFrames.DataFrame
        @test size(tmp) == (237,3)
        tmp = heterozygosity(cats, by = "population")
        @test tmp isa DataFrames.DataFrame
        @test size(tmp) == (17,4)
        @test samplehet(cats, "N215") == 1/3
        @test_throws ArgumentError samplehet(cats, "M115")
    end
end

end #module