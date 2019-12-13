module PopGen

##   o O       o O       o O       o O       o O
## o | | O   o | | O   o | | O   o | | O   o | | O
## | | | | O | | | | Dependencies| | | | O | | | | O
## O | | o   O | | o   O | | o   O | | o   O | | o
##   O o       O o       O o       O o       O o


using DataFrames,
      PlotlyJS,
      GeneticVariation,
      Distributions,
      MultipleTesting,
      LinearAlgebra,
      Convex,
      ECOS,
      ProgressMeter


export PopObj,
    summary,
    nancycats,
    gulfsharks,
    csv,
    genepop,
    vcf,
    samples,
    loci,
    isolate_genotypes,
    locations, locations!,
    population, populations, population!, populations!,
    remove_inds!,
    remove_loci!,
    missing,
    heterozygosity, het, He,
    hwe_test, hwe,
    plot_missing,
    plot_locations



##   o O       o O       o O       o O       o O
## o | | O   o | | O   o | | O   o | | O   o | | O
## | | | | O | | | |Include Files| | | | O | | | | O
## O | | o   O | | o   O | | o   O | | o   O | | o
##   O o       O o       O o       O o       O o

include("PopObj.jl")
include("Read.jl")
include("Manipulate.jl")
include("Plotting.jl")
include("Datasets.jl")
include("AlleleFreq.jl")
include("HardyWeinberg.jl")

end # module PopGen
