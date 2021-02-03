module  TestAlleleFreqs

using PopGen
using DataFrames
using Test

cats = @nancycats;
sharks = @gulfsharks;

@testset "PopGen.allele_freq" begin
    @test length(PopGen.allele_freq(cats)) == 9
    @test eltype(PopGen.allele_freq(cats)) == Dict{Int16,Float64}
    @test typeof(PopGen.allele_freq(cats.loci.genotype)) == Dict{Int16,Float64}
    df_cats = DataFrames.combine(
        groupby(cats.loci, :population),
        :genotype => PopGen.allele_freq => :genos
    )
    @test size(df_cats) == (17,2)
    @test eltype(df_cats.genos) == Dict{Int16,Float64}
    @test length(PopGen.allele_freq(cats, "fca8")) == 16
    @test typeof(PopGen.allele_freq(cats, "fca8")) == Dict{Int16,Float64}
    
    @test length(PopGen.allele_freq(sharks)) == 2209
    @test eltype(PopGen.allele_freq(sharks)) == Dict{Int8,Float64}
    @test typeof(PopGen.allele_freq(sharks.loci.genotype)) == Dict{Int8,Float64}
    df_sharks = DataFrames.combine(
        groupby(sharks.loci, :population),
        :genotype => PopGen.allele_freq => :genos
    )
    @test size(df_sharks) == (7,2)
    @test eltype(df_sharks.genos) == Dict{Int8,Float64}
    @test length(PopGen.allele_freq(sharks, "contig_2784")) == 2
    @test typeof(PopGen.allele_freq(sharks, "contig_2784")) == Dict{Int8,Float64}
end

@testset "other PopGen.allele_freqs" begin
    @test PopGen.allele_freq_vec([(1,1), (2,2), (2,1)]) == [0.5, 0.5]
    @test PopGen.allele_freq_vec(missing) === missing
end

@testset "geno counts" begin
    @test length(PopGen.geno_count_observed(cats.loci.genotype)) == 295
    @test typeof(PopGen.geno_count_observed(cats.loci.genotype)) <: Dict{<:Tuple,Int64}
    @test length(PopGen.geno_count_expected(cats.loci.genotype)) == 6241
    @test typeof(PopGen.geno_count_expected(cats.loci.genotype)) <: Dict{<:Tuple,Float64}

    @test length(PopGen.geno_count_observed(sharks.loci.genotype)) == 12
    @test typeof(PopGen.geno_count_observed(sharks.loci.genotype)) <: Dict{<:Tuple,Int64}
    @test length(PopGen.geno_count_expected(sharks.loci.genotype)) == 25
    @test typeof(PopGen.geno_count_expected(sharks.loci.genotype)) <: Dict{<:Tuple,Float64}
end

@testset "geno freqs" begin
    @test length(PopGen.geno_freq(cats.loci.genotype)) == 295
    @test typeof(PopGen.geno_freq(cats.loci.genotype)) <: Dict{<:Tuple,Float64}
    @test length(PopGen.geno_freq(cats, "fca8")) == 51
    @test typeof(PopGen.geno_freq(cats, "fca8")) <: Dict{<:Tuple,Float64}
    @test typeof(PopGen.geno_freq(cats, "fca8", population = true)) == DataFrame
    @test size(PopGen.geno_freq(cats, "fca8", population = true)) == (17,2)

    @test length(PopGen.geno_freq_expected(cats.loci.genotype)) == 6241
    @test typeof(PopGen.geno_freq_expected(cats.loci.genotype)) <: Dict{<:Tuple,Float64}
    @test length(PopGen.geno_freq_expected(cats, "fca8")) == 256
    @test typeof(PopGen.geno_freq_expected(cats, "fca8")) <: Dict{<:Tuple,Float64}
    @test typeof(PopGen.geno_freq_expected(cats, "fca8", population = true)) == DataFrame
    @test size(PopGen.geno_freq_expected(cats, "fca8", population = true)) == (17,2)

    @test length(PopGen.geno_freq(sharks.loci.genotype)) == 12
    @test length(PopGen.geno_freq(sharks, "contig_2784")) == 2
    @test typeof(PopGen.geno_freq(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
    @test typeof(PopGen.geno_freq(sharks, "contig_2784", population = true)) == DataFrame
    @test size(PopGen.geno_freq(sharks, "contig_2784", population = true)) == (7,2)

    @test length(PopGen.geno_freq_expected(sharks.loci.genotype)) == 25
    @test typeof(PopGen.geno_freq_expected(sharks.loci.genotype)) <: Dict{<:Tuple,Float64}
    @test length(PopGen.geno_freq_expected(sharks, "contig_2784")) == 4
    @test typeof(PopGen.geno_freq_expected(sharks, "contig_2784")) <: Dict{<:Tuple,Float64}
    @test typeof(PopGen.geno_freq_expected(sharks, "contig_2784", population = true)) == DataFrame
    @test size(PopGen.geno_freq_expected(sharks, "contig_2784", population = true)) == (7,2)
end

end # module