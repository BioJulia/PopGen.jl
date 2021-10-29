"""
# Population genetics analyses in Julia
Repository:    https://www.github.com/biojulia/PopGen.jl/
Documentation: https://biojulia.net/PopGen.jl/

\nA few things things you can do to get started:

## Import Data
- `PopGen.read(filename; kwargs...)`
- `genepop(infile; kwargs...)`  or similar file-specific importer
- use available `@gulfsharks` or `@nancycats` datasets

## Explore PopData
- `populations(PopData)` to view population information
- `loci(PopData)` to view locus names
- `samplenames(PopData)` to view sample names
- `missingdata(PopData, by = ...)` to view missing information

## Manipulate PopData
- `populations!(PopData, ...)` to rename populations
- `locations!(PopData, ...)` to add geographical coordinates
- `exclude!(PopData; kwargs...)` to selectively remove data

## Analyses
- `richness(PopData)` to calculate allelic richness
- `relatedness(PopData, method = ...)` to get pairwise relatedness of individuals
- `summary(PopData)` to calculate F-statistics, heterozygosity, etc.
- `hwetest(PopData)` to test for Hardy-Weinberg Equilibrium
- `pairwisefst(PopData)` to calculate FST between pairs of populations
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
    export genodata, metadata, info
    export isbiallelic, ishom, ishet
    export delimited, csv, genepop, vcf, bcf, @nancycats, @gulfsharks
    export sampleinfo, sampleinfo!, locusinfo, locusinfo!, samplenames, loci
    export copy, size, sort, dropmonomorphic, dropmonomorphic!
    export dropmultiallelic, dropmultiallelic!
    export locationdata, locationdata!
    export locidataframe, locimatrix
    export genotypes, genotypes
    export populations, populations!
    export exclude, remove, omit, exclude!, remove!, omit!, keep, keep!
    export filter, filter!
end
@reexport import PopGenCore: read, write

using Distributions, DataFrames, PooledArrays
using Random: shuffle
using ProgressMeter
using MultipleTesting, StatsBase


#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | |Include Files| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

include("Utils.jl")
# heterozygosity functions
include("Heterozygosity.jl")
# manipulation and exploration
include("DataExploration.jl")
# summary information
include("SummaryInfo.jl")
#Analyses
include("HardyWeinberg.jl")
include("FStatistics/FstGlobal.jl")
include("FStatistics/FstByLocus.jl")
include("FStatistics/PairwiseFST.jl")
include("FStatistics/FstPermutations.jl")
include("Relatedness/PairwiseRelatedness.jl")
include("Relatedness/RelatednessMoments.jl")
include("Relatedness/RelatednessPostHocs.jl")

end # module PopGen
