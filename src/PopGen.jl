module PopGen

##   o O       o O       o O       o O       o O
## o | | O   o | | O   o | | O   o | | O   o | | O
## | | | | O | | | | Dependencies| | | | O | | | | O
## O | | o   O | | o   O | | o   O | | o   O | | o
##   O o       O o       O o       O o       O o

using CSV, JuliaDB, JuliaDBMeta, MultipleTesting, Random, StatsBase, CategoricalArrays

#=
using Convex,
      CSV,
      Distributions,
      ECOS,
      GeneticVariation,
      JuliaDB,
      JuliaDBMeta,
      LinearAlgebra,
      MultipleTesting,
      ProgressMeter,
      Random,
      StatsBase
=#

export PopObj,
    summary,
    nancycats,
    gulfsharks,
    delimited, csv,
    genepop,
    bcf, vcf,
    samples,
    loci,
    isolate_genotypes,
    locations, locations!,
    population, populations, population!, populations!, popnames!,
    relatedness, pairwise_relatedness, kinship,
    remove_inds!,
    remove_loci!,
    missing,
    show,
    heterozygosity, het, He,
    hwe_test, hwe
    #plot_missing,
    #plot_locations



##   o O       o O       o O       o O       o O
## o | | O   o | | O   o | | O   o | | O   o | | O
## | | | | O | | | |Include Files| | | | O | | | | O
## O | | o   O | | o   O | | o   O | | o   O | | o
##   O o       O o       O o       O o       O o

# The types
include("Types.jl")
# file io
include("io/ioUtils.jl")
include("io/Delimited.jl")
include("io/Genepop.jl")
include("io/VariantCall.jl")
# example data
include("Datasets.jl")
# manipulation commands
include("Manipulate.jl")
# allele frequency and heterozygosity functions
include("AlleleFreq.jl")
include("HardyWeinberg.jl")
# summary information
include("SummaryInfo.jl")
#include("PairwiseRelatedness.jl")  # not yet ready
#include("PlotRecipes.jl")  # not yet ready
end # module PopGen
