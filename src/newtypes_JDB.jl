using DataFrames, BenchmarkTools, JuliaDBMeta, DataFrames, PopGen, CategoricalArrays, DataFramesMeta, StatsBase
import JuliaDB
#=
abstract type PopObj end

mutable struct PopSample <: PopObj
    name::String
    population::String
    ploidy::Int8
    longitude::Union{Missing, Float32}
    latitude::Union{Missing, Float32}
end

struct PopData <: PopObj
    samples::Vector{PopSample}
    loci::T where T<:IndexedTable
end
=#

##### gulfsharks stuff
function gulfsharks_lf_table()
  sharks = gulfsharks();
  sharks_df =  insertcols!(sharks.loci,1,:name => sharks.samples.name)
  sharks_df =  insertcols!(sharks.loci,2,:population => sharks.samples.population)
  sharks_df_long = DataFrames.stack(sharks_df, 3:2215)
  rename!(sharks_df_long, [:locus, :genotype, :name, :population])
  sharks_df_long.locus = string.(sharks_df_long.locus)
  sharks_table_noncat = sharks_df_long |> table
  #=Base.summarysize(sharks_df_long)
  12671248 =#
  categorical!(sharks_df_long, :locus)
  categorical!(sharks_df_long, :population)
  categorical!(sharks_df_long, :name)


  sharks_table = sharks_df_long |> table
  return sharks_df, sharks_df_long, sharks_table, sharks_table_noncat
end
a,b,c,d = gulfsharks_lf_table()
#a standard loci DF
#b long-format loci DF
#c long-format JDB table with CategoricalArrays
#d long-format JDB table without categorical

sharks = gulfsharks();

### missing by pop
@btime @groupby c :population {missing = sum(ismissing.(:genotype))}
@btime @groupby d :population {missing = sum(ismissing.(:genotype))}
@btime PopGen.missing(sharks, mode = "pop")
@btime by(b, :population, miss = :genotype => i -> sum(ismissing.(i)))

### missing by loc
@btime @groupby c :locus {missing = sum(ismissing.(:genotype))}
@btime @groupby d :locus {missing = sum(ismissing.(:genotype))}
@btime PopGen.missing(sharks, mode = "locus")
@btime @by(b, :locus, miss = sum(ismissing.(:genotype)))
@btime by(b, :locus, miss = :genotype => i -> sum(ismissing.(i)))



# het?
@btime by(b, :population, het = :genotype => i -> mean(length.(unique.(skipmissing(i) |> collect))))
# 156.256 ms (2787474 allocations: 236.78 MiB)
@btime by(b, :locus, het = :genotype => i -> mean(length.(unique.(skipmissing(i) |> collect))))
# 179.847 ms (3055436 allocations: 244.91 MiB)

@btime @groupby c (:locus) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 543.717 ms (2851555 allocations: 241.94 MiB)
@btime @groupby c (:locus, :population) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 2.871 s (3184984 allocations: 306.68 MiB)
@btime @groupby d (:locus) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 246.636 ms (2844834 allocations: 241.48 MiB)
@btime @groupby d (:population) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 163.241 ms (2787463 allocations: 230.18 MiB)
@btime @groupby d (:locus, :population) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 278.037 ms (3149464 allocations: 306.51 MiB)

@btime PopGen.het_observed(sharks)
# 1.318 s (1894933 allocations: 48.14 MiB)
@btime PopGen.het_population_obs(sharks)
# 71.318 ms (580830 allocations: 28.91 MiB)


# DataFrames
## CategoricalArrays
@btime by(b, :locus, het = :genotype => i -> mean(length.(unique.(skipmissing(i) |> collect))))
# 179.847 ms (3055436 allocations: 244.91 MiB)

## Regular String arrays
@btime by(a, :locus, het = :genotype => i -> mean(length.(unique.(skipmissing(i) |> collect))))
#  175.727 ms (2830513 allocations: 245.18 MiB)

#JuliaDB
## CategoricalArrays
@btime @groupby $c (:locus) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 543.717 ms (2851555 allocations: 241.94 MiB)

## Regular String arrays
@btime @groupby $d (:locus) {hetero = mean(length.(unique.(skipmissing(:genotype) |> collect)) .> 1)}
# 246.636 ms (2844834 allocations: 241.48 MiB)
