"""
# Quickstart for PopGen
Documentation: https://pdimens.github.io/PopGen.jl/

Motivational(?) quote: $(motivational_quote())

\nA few things things you can do to get started:

## Load in data
- `read_from(filename; kwargs...)`
- `genepop(infile; kwargs...)`  or similar file-specific importer
- use available `gulfsharks()` or `nancycats()` datasets

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

using CSV, Distributions, DataFrames
using FileIO, JLD2
using MultipleTesting, Random, StatsBase
using GeneticVariation: VCF, BCF, header

#=
using Convex,
      ECOS,
      LinearAlgebra,
      ProgressMeter,
=#

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
include("io/VariantCall.jl")
# example data
include("io/Datasets.jl")
# utility functions
include("Utils.jl")
include("Permutations.jl")
# manipulation and exploration
include("DataExploration.jl")
include("Manipulate.jl")
# allele frequency and heterozygosity functions
include("AlleleFreq.jl")
include("Heterozygosity.jl")
# summary information
include("SummaryInfo.jl")
#Analyses
include("HardyWeinberg.jl")
#include(PairwiseRelatedness.jl)  # not yet ready
#include(PlotRecipes.jl)  # not yet ready

end # module PopGen
