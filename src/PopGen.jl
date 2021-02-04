"""
# Quickstart for PopGen
Documentation: https://pdimens.github.io/PopGen.jl/

Motivational(?) quote: $(motivational_quote())

\nA few things things you can do to get started:

## Load in data
- `read_from(filename; kwargs...)`
- `genepop(infile; kwargs...)`  or similar file-specific importer
- use available `@gulfsharks` or `@nancycats` datasets

## Explore PopData
- `populations(PopData)` to view population information
- `loci(PopData)` to view locus names
- `samples(PopData)` to view sample names
- `missing(PopData, by = ...)` to view missing information

## Manipulate PopData
- `populations!(PopData, ...)` to rename populations
- `locations!(PopData, ...)` to add geographical coordinates
- `exclude!(PopData; kwargs...)` to selectively remove data

## Analyses
- `richness(PopData)` to calculate allelic richness
- `allele_avg(PopData)` to calculate average # of alleles
- `summary(PopData)` to calculate F-statistics, heterozygosity, etc.
- `hwe_test(PopData)` to test for Hardy-Weinberg Equilibrium
"""
module PopGen

#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | | Dependencies| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

using CSV, Distributions, DataFrames, PooledArrays, StaticArrays
using FileIO, JLD2, Requires, ProgressMeter
using MultipleTesting, Random, StatsBase
using RandomNumbers.Xorshifts


#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | |Include Files| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

# the types
include("Types.jl")
# file io
include("io/ioUtils.jl")
include("io/Delimited.jl")
include("io/Genepop.jl")
include("io/Read.jl")
include("io/Structure.jl")
@init @require GeneticVariation="9bc6ac9d-e6b2-5f70-b0a8-242a01662520" begin
    include("io/VariantCall.jl")
end
@init @require GeneticVariation="9bc6ac9d-e6b2-5f70-b0a8-242a01662520" begin
    @require GZip="92fee26a-97fe-5a0c-ad85-20a5f3185b63" include("io/VariantCallGz.jl")
end
# example data
include("io/Datasets.jl")
# utility functions
include("Utils.jl")
include("Permutations.jl")
# allele frequency and heterozygosity functions
include("AlleleFreq.jl")
include("Heterozygosity.jl")
# manipulation and exploration
include("DataExploration.jl")
include("Manipulate.jl")
# summary information
include("SummaryInfo.jl")
#Analyses
include("HardyWeinberg.jl")
include("Relatedness/PairwiseRelatedness.jl")
include("Relatedness/RelatednessMoments.jl")
include("Relatedness/RelatednessPostHocs.jl")
#include(PlotRecipes.jl)  # not yet ready

end # module PopGen
