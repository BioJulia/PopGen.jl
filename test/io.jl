module  TestIO

using PopGen
using GeneticVariation, GZip
using Test

cats_gen = normpath(joinpath(@__DIR__,"..","data/source", "nancycats.gen"))
sharks_gen = normpath(joinpath(@__DIR__,"..","data/source", "gulfsharks.gen"))
oyster_vcf = normpath(joinpath(@__DIR__,"..","data/source", "filtered_oyster.vcf"))

@testset "Genepop io" begin
    @test typeof(read_from(cats_gen, silent = true)) == PopData
    @test typeof(read_from(sharks_gen, silent = true)) == PopData
    @test typeof(genepop(cats_gen, silent = true)) == PopData
    @test typeof(genepop(sharks_gen, silent = true)) == PopData
end


@testset "VCF io" begin
    @test typeof(read_from(oyster_vcf, silent = true)) == PopData
    @test typeof(vcf(oyster_vcf, silent = true)) == PopData    
end

end # module

