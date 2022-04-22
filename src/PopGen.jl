"""
# Population genetics analyses in Julia
Repository:    https://www.github.com/biojulia/PopGen.jl/
Documentation: https://biojulia.net/PopGen.jl/

\nA few things things you can do to get started:

## Import Data
- `PopGen.read(filename, kwargs...)`
- `genepop(infile, kwargs...)`  or similar file-specific importer
- use available `@gulfsharks` or `@nancycats` datasets

## Explore PopData
- `populations(PopData)` to view population information
- `loci(PopData)` to view locus names
- `samplenames(PopData)` to view sample names
- `missingdata(PopData, by = ...)` to view missing information

## Manipulate PopData
- `populations!(PopData, ...)` to rename populations
- `locations!(PopData, ...)` to add geographical coordinates
- `exclude!(PopData, kwargs...)` to selectively remove data

## Analyses
- `richness(PopData)` to calculate allelic richness
- `Kinship(PopData, method = ...)` to get pairwise Kinship of individuals
- `summary(PopData)` to calculate F-statistics, heterozygosity, etc.
- `hwetest(PopData)` to test for Hardy-Weinberg Equilibrium
- `pairwisefst(PopData)` to calculate FST between pairs of populations
- `kmeans(PopData, 2:k)` to perform Kmeans++ clustering
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

using Distributions, DataFrames, PooledArrays NamedArrays
using Distances
using Random: shuffle
using ProgressMeter
using MultipleTesting, StatsBase
import Clustering: kmeans, kmedoids, hclust, Hclust, cutree, fuzzy_cmeans, dbscan
import MultivariateStats: fit, PCA


#   o O       o O       o O       o O       o O
# o | | O   o | | O   o | | O   o | | O   o | | O
# | | | | O | | | |Include Files| | | | O | | | | O
# O | | o   O | | o   O | | o   O | | o   O | | o
#   O o       O o       O o       O o       O o

include("Utils.jl")
include("AlleleMatrices.jl")

# heterozygosity functions
include("Heterozygosity.jl")
export heterozygosity, samplehet

# manipulation and exploration
include("DataExploration.jl")
export pairwiseidentical, missingdata, genofreqtable, allelefreqtable

# summary information
include("SummaryInfo.jl")
export alleleaverage, richness, summary, summarystats

#Analyses
include("HardyWeinberg.jl")
export hwetest, hwe

include("FStatistics/FstGlobal.jl")
include("FStatistics/FstByLocus.jl")
include("FStatistics/PairwiseFST.jl")
export pairwisefst

include("FStatistics/FstPermutations.jl")
include("Kinship/KinshipPairwise.jl")
export kinship, mergeKinship

include("Kinship/KinshipMoments.jl")
export QuellerGoodnight, Ritland, Lynch, LynchRitland, LynchLi, LiHorvitz, Moran, Blouin, Loiselle #, Wang

include("Kinship/KinshipPostHocs.jl")
export kinshipposthoc

include("Clustering.jl")
export cluster, kmeans, kmedoids, hclust, cutree, fuzzycmeans, dbscan

include("PCA.jl")
export pca

end # module PopGen
