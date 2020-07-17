
export pairwise_identical, missing

"""
    missing(data::PopData; by::String = "sample")
Get missing genotype information in a `PopData`. Specify a mode of operation
to return a DataFrame corresponding with that missing information.

#### Modes
- "sample" - returns a count and list of missing loci per individual (default)
- "pop" - returns a count of missing genotypes per population
- "locus" - returns a count of missing genotypes per locus
- "full" - returns a count of missing genotypes per locus per population

### Example:
```
missing(gulfsharks(), by = "pop")
```
"""

@inline function Base.missing(data::PopData; by::String = "sample")
    if by ∈ ["sample", "individual"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, :name),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    elseif by ∈ ["pop", "population"]
                DataFrames.combine(
            DataFrames.groupby(data.loci, :population),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    elseif by ∈ ["locus", "loci"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, :locus),
            :genotype => (i -> count(ismissing, i)) => :missing
        )

    elseif by ∈ ["detailed", "full"]
        DataFrames.combine(
            DataFrames.groupby(data.loci, [:locus, :population]),
            :genotype => (i -> count(ismissing, i)) => :missing
        )
    else
        @error "Mode \"$by\" not recognized. Please specify one of: sample, pop, locus, or full"
        missing(data)
    end
end

#TODO add to docs (Data Exploration page and API)
function pairwise_identical(data::PopData)
    sample_names = samples(data)
    sample_pairs = [tuple(sample_names[i], sample_names[j]) for i in 1:length(sample_names)-1 for j in i+1:length(sample_names)]
    n = length(sample_pairs)
    perc_ident_vec = Vector{Float64}(undef, n)
    n_vec = Vector{Int}(undef, n)
    idx = 0
    @inbounds for (sample_n, sample_1) in enumerate(sample_names[1:end-1])
        geno_1 = get_genotypes(data, sample_1)
        len_1 = length(collect(skipmissing(geno_1)))
        @inbounds Base.Threads.@threads for sample_2 in sample_names[sample_n+1:end]
            idx += 1
            geno_2 = get_genotypes(data, sample_2)
            len_2 = length(collect(skipmissing(geno_2)))
            shared_geno = minimum([len_1, len_2])
            shared = sum(skipmissing(geno_1 .== geno_2))
            perc_ident_vec[idx] = round(shared/shared_geno, digits = 2)
            n_vec[idx] = shared_geno
        end
    end
    DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :identical => perc_ident_vec, :n => n_vec)
end