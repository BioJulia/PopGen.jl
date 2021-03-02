module  TestIO

using PopGen
using GeneticVariation, GZip
using Test

cats_gen = normpath(joinpath(@__DIR__,"..","data/", "nancycats.gen"))
sharks_csv = normpath(joinpath(@__DIR__,"..","data/", "gulfsharks.csv"))
oyster_vcf = normpath(joinpath(@__DIR__,"..","data/", "filtered_oyster.vcf"))

@testset "Genepop io" begin
    @test typeof(read_from(cats_gen, silent = true)) == PopData
    @test typeof(genepop(cats_gen, silent = true)) == PopData
end

@testset "delimited io" begin
    @test typeof(read_from(sharks_csv, silent = true)) == PopData
    @test typeof(delimited(sharks_csv, silent = true)) == PopData
end

@testset "VCF io" begin
    @test typeof(read_from(oyster_vcf, silent = true)) == PopData
    @test typeof(vcf(oyster_vcf, silent = true)) == PopData    
end

end # module

