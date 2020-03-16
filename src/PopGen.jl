module PopGen

#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | | Dependencies| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

using CategoricalArrays, CSV, Distributions, JuliaDB, JuliaDBMeta, MultipleTesting, Random, StatsBase
#import IndexedTables: reindex

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
    exclude_loci, omit_loci, remove_loci,
    exclude_samples, omit_samples, remove_samples,
    genepop,
    bcf, vcf,
    samples,
    loci,
    locus,
    isolate_genotypes,
    locations, locations!,
    population, populations, population!, populations!, popnames!,
    relatedness, pairwise_relatedness, kinship,
    meta,
    missing,
    reindex,
    show,
    heterozygosity, het, He,
    hwe_test, hwe
    #plot_missing,
    #plot_locations



#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | |Include Files| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

# The types
include("Types.jl")
# file io
include("io/ioUtils.jl")
#include("io/Delimited.jl")
include("io/Genepop.jl")
#include("io/VariantCall.jl")
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
