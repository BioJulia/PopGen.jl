module  TestTypesDatasets

using PopGen
using DataFrames
using PooledArrays
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "Dataset types" begin
    @test typeof(cats) == PopData
    @test typeof(cats.meta) == DataFrame
    @test typeof(cats.loci) == DataFrame
    @test typeof(sharks) == PopData
    @test typeof(sharks.meta) == DataFrame
    @test typeof(sharks.loci) == DataFrame
end

@testset "Dataset dimensions" begin
    @test size(cats.meta) == (237,5)
    @test size(cats.loci) == (2133,4)
    @test size(sharks.meta) == (212,5)
    @test size(sharks.loci) == (468308,4)
end

@testset "Dataset column names" begin
    meta_names_sorted = ["latitude","longitude","name","ploidy","population"]
    loci_names_sorted = ["genotype","locus","name","population"]
    @test sort(names(cats.meta)) == meta_names_sorted
    @test sort(names(cats.loci)) == loci_names_sorted
    @test sort(names(sharks.meta)) == meta_names_sorted
    @test sort(names(sharks.loci)) == loci_names_sorted
end

@testset "Nancycats column types" begin
    @test typeof(cats.meta.name) == Vector{String}
    @test typeof(cats.meta.population) == Vector{String}
    @test typeof(cats.meta.ploidy) == Vector{Int8}
    @test typeof(cats.meta.latitude) ==  Vector{Union{Missing, Float32}}
    @test typeof(cats.meta.longitude) == Vector{Union{Missing, Float32}}
    @test typeof(cats.loci.name) <: PooledArray
    @test eltype(cats.loci.name) == String
    @test typeof(cats.loci.population) <: PooledArray
    @test eltype(cats.loci.population) == String
    @test typeof(cats.loci.locus) <: PooledArray
    @test eltype(cats.loci.locus) == String
    @test typeof(cats.loci.genotype) <: GenoArray
    @test eltype(cats.loci.genotype) <: Union{Missing, Genotype}
end

@testset "Gulfsharks column types" begin
    @test typeof(sharks.meta.name) == Vector{String}
    @test typeof(sharks.meta.population) == Vector{String}
    @test typeof(sharks.meta.ploidy) == Vector{Int8}
    @test typeof(sharks.meta.latitude) ==  Vector{Float64}
    @test typeof(sharks.meta.longitude) == Vector{Float64}
    @test typeof(sharks.loci.name) <: PooledArray
    @test eltype(sharks.loci.name) == String
    @test typeof(sharks.loci.population) <: PooledArray
    @test eltype(sharks.loci.population) == String
    @test typeof(sharks.loci.locus) <: PooledArray
    @test eltype(sharks.loci.locus) == String
    @test typeof(sharks.loci.genotype) <: GenoArray
    @test eltype(sharks.loci.genotype) <: Union{Missing, Genotype}
end


end  # module