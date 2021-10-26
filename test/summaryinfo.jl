module  SummaryInfo

using PopGen
using DataFrames
using Test

cats = @nancycats;

@testset "alleleaverage" begin
    @test length(alleleaverage(cats,  rounding = true)) == 2
    @test alleleaverage(cats,  rounding = true) isa NamedTuple
    @test length(alleleaverage(cats,  rounding = false)) == 2
    @test alleleaverage(cats,  rounding = false) isa NamedTuple
end

@testset "richness" begin
    @test richness(cats, by = "locus") isa DataFrame
    @test size(richness(cats, by = "locus")) == (9, 2)
    @test richness(cats, by = "population") isa DataFrame
    @test size(richness(cats, by = "population")) == (153, 3)
end

@testset "summary F/D/etc. stats" begin
    smry_glob = summary(cats)
    smry_loc = summary(cats, by = "locus")
    @test smry_glob isa DataFrame
    @test smry_loc isa DataFrame
    @test size(smry_glob) == (1, 10)
    @test size(smry_loc) ==(9, 11)
end

end # module