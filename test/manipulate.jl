module  TestManipulate

using PopGen
using DataFrames
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "locations" begin
    x = rand(length(samples(cats))) ; y = rand(length(samples(cats)))
    locations!(cats, long = x, lat = y)    
    @test cats.meta.longitude == Float32.(x)
    @test cats.meta.latitude == Float32.(y)
    @test locations(cats).longitude == Float32.(x)
    @test locations(cats).latitude == Float32.(y)
end

@testset "decimal-minutes locations" begin
    x = fill("11 22.33W", length(samples(cats))) ; y = fill("-41 31.52", length(samples(cats)))
    locations!(cats, long = x, lat = y)
    @test all(cats.meta.longitude .== Float32(-11.3722))
    @test all(cats.meta.latitude .== Float32(-41.5253))
end

@testset "loci and genotypes" begin
    @test length(loci(cats)) == 9
    @test get_genotype(cats, sample = "N115", locus = "fca8") == (135,135)
    N115 = get_genotypes(cats, "N115")
    @test length(N115) == 9
    @test typeof(N115) == Vector{Union{Missing, Tuple{Int16,Int16}}}
    @test typeof(get_genotypes(cats, sample = "N115" , locus = "fca8")) <: SubDataFrame
    @test names(get_genotypes(cats, sample = "N115" , locus = "fca8")) == ["name", "population", "locus", "genotype"]
    @test size(get_genotypes(cats, sample = ["N115", "N7"] , locus = "fca8")) == (2,4)
    @test size(get_genotypes(cats, sample = "N115" , locus = ["fca8", "fca37"])) == (2,4)
    @test size(get_genotypes(cats, sample = ["N115", "N7"] , locus = ["fca8", "fca37"])) == (4,4)
    @test length(genotypes(sharks, "contig_475")) == 212
end

@testset "populations" begin
    @test size(populations(cats)) == (17,2)
    @test typeof(populations(cats).population) == Vector{String}
    @test typeof(populations(cats).count) == Vector{Int}

    rn_dict = Dict("1" => "one", "2" => "two")
    populations!(cats, rn_dict)
    @test "one" ∈ cats.meta.population && "1" ∉ cats.meta.population
    @test "two" ∈ cats.meta.population && "2" ∉ cats.meta.population

    rn_vect = string.(1:17)
    populations!(cats, rn_vect)
    @test "one" ∉ cats.meta.population && "1" ∈ cats.meta.population
    @test "two" ∉ cats.meta.population && "2" ∈ cats.meta.population
 
    rn_vect_ii = ["N215", "N297"]
    rn_vect_iin = ["one", "seventeen"]
    populations!(cats, rn_vect_ii, rn_vect_iin)
    @test "one" ∈ cats.meta.population && "seventeen" ∈ cats.meta.population
    @test cats.meta[cats.meta.name .== "N215", :population] == ["one"]
    @test cats.meta[cats.meta.name .== "N297", :population] == ["seventeen"]
end

@testset "exclusion" begin
    tmp = exclude(nancycats(), name = "N100", population = ["1", "15"])
    @test length(samples(tmp)) == 215
    @test size(populations(tmp)) == (15,2)

    tmp = exclude(nancycats(), names = "N102", loci = "fca8", population = "3")
    @test length(loci(tmp)) == 8
    @test size(populations(tmp)) == (16,2)
    @test length(samples(tmp)) == 225
end

end # module