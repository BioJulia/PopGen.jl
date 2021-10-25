module  TestDataExploration

using PopGen
using DataFrames
using Test

cats = @nancycats;

@testset "Data Exploration" begin

    @testset "missing data" begin
        @test typeof(missingdata(cats)) == DataFrame
        @test typeof(missingdata(cats, by = "population")) == DataFrame
        @test typeof(missingdata(cats, by = "locus")) == DataFrame
        @test typeof(missingdata(cats, by = "full")) == DataFrame
    end

    @testset "pairwise identical" begin
        pw_i = pairwiseidentical(cats)
        @test typeof(pw_i) == DataFrame
        @test size(pw_i) == (27966, 4)
        pw_i_2 = pairwiseidentical(cats, cats.sampleinfo.name[1:10])
        @test typeof(pw_i_2) == DataFrame
        @test size(pw_i_2) == (45, 4)
    end

    @testset "frequency tables" begin
        @test typeof(genofreqtable(cats)) == DataFrame
        @test typeof(genofreqtable(cats, by = "population")) == DataFrame
        @test typeof(allelefreqtable(cats)) == DataFrame
        @test typeof(allelefreqtable(cats, by = "population")) == DataFrame
    end
end

end # module