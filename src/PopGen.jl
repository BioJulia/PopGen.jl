module PopGen

#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | | Dependencies| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

using CategoricalArrays, CSV, Distributions, IndexedTables, JuliaDBMeta, LoopVectorization, MultipleTesting, Random, StatsBase

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
#include("io/VariantCall.jl")
# utility functions
include("Utils.jl")
# example data
include("Datasets.jl")
# manipulation commands
include("Manipulate.jl")
# allele frequency and heterozygosity functions
include("AlleleFreq.jl")
include("Heterozygosity.jl")
# summary information
include("SummaryInfo.jl")
#Analyses
include("HardyWeinberg.jl")
#include("PairwiseRelatedness.jl")  # not yet ready
#include("PlotRecipes.jl")  # not yet ready

end # module PopGen
