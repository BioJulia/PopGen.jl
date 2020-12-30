module  TestDataExploration

using PopGen
using DataFrames
using Test

cats = @nancycats;

@testset "missing()" begin
    @test typeof(missing_data(cats)) == DataFrame
    @test typeof(missing_data(cats, by = "population")) == DataFrame
    @test typeof(missing_data(cats, by = "locus")) == DataFrame
    @test typeof(missing_data(cats, by = "full")) == DataFrame
end

@testset "pairwise_identical" begin
    pw_i = pairwise_identical(cats)
    @test typeof(pw_i) == DataFrame
    @test size(pw_i) == (27966, 4)
    pw_i_2 = pairwise_identical(cats, samples(cats)[1:10])
    @test typeof(pw_i_2) == DataFrame
    @test size(pw_i_2) == (45, 4)
end

@testset "frequency tables" begin
    @test typeof(geno_freqtable(cats)) == DataFrame
    @test typeof(geno_freqtable(cats, by = "population")) == DataFrame
    @test typeof(allele_freqtable(cats)) == DataFrame
    @test typeof(allele_freqtable(cats, by = "population")) == DataFrame
end

end # module