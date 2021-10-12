"""
# Population genetics analyses in Julia
Documentation: https://pdimens.github.io/PopGen.jl/

\nA few things things you can do to get started:

## Load in data
- `PopGen.read(filename; kwargs...)`
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
- `relatedness(PopData, method = ...)` to get pairwise relatedness of individuals
- `summary(PopData)` to calculate F-statistics, heterozygosity, etc.
- `hwe_test(PopData)` to test for Hardy-Weinberg Equilibrium
"""
module PopGen

#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | | Dependencies| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

using Reexport
using PopGenCore
@reexport module PopGenCore
    export PopData, PopDataInfo, GenoArray, Genotype, SNP, Msat
    export isbiallelic, ishom, ishet
    export delimited, csv, genepop, vcf, bcf, @nancycats, @gulfsharks
    export ishom, ishet
end


using Distributions, DataFrames, PooledArrays
using ProgressMeter
using MultipleTesting, StatsBase


#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | |Include Files| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

# allele frequency and heterozygosity functions
include("AlleleFreq.jl")
include("Heterozygosity.jl")
# manipulation and exploration
include("DataExploration.jl")
# summary information
include("SummaryInfo.jl")
#Analyses
include("HardyWeinberg.jl")
include("FStatistics/FstMethods.jl")
include("FStatistics/PairwiseFST.jl")
include("FStatistics/FstPermutations.jl")
include("Relatedness/PairwiseRelatedness.jl")
include("Relatedness/RelatednessMoments.jl")
include("Relatedness/RelatednessPostHocs.jl")

end # module PopGen
