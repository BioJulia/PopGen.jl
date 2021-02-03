module  SummaryInfo

using PopGen
using DataFrames
using Test

cats = @nancycats;

@testset "allele_avg" begin
    @test length(allele_avg(cats,  rounding = true)) == 2
    @test typeof(allele_avg(cats,  rounding = true)) <: NamedTuple
    @test length(allele_avg(cats,  rounding = false)) == 2
    @test typeof(allele_avg(cats,  rounding = false)) <: NamedTuple
end

@testset "richness" begin
    @test typeof(richness(cats, by = "locus")) == DataFrame
    @test size(richness(cats, by = "locus")) == (9, 2)
    @test typeof(richness(cats, by = "population")) == DataFrame
    @test size(richness(cats, by = "population")) == (153, 3)
end

@testset "summary F/D/etc. stats" begin
    smry_glob = summary(cats)
    smry_loc = summary(cats, by = "locus")
    @test typeof(smry_glob) == DataFrame
    @test typeof(smry_loc) == DataFrame
    @test size(smry_glob) == (1, 10)
    @test size(smry_loc) ==(9, 11)
end

end # module