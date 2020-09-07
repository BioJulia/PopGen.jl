
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


"""
    pairwise_identical(data::PopData)
Return a table of the percent of identical genotypes at each locus between pairs of individuals.
"""
function pairwise_identical(data::PopData)
    sample_names = collect(samples(data))
    pairwise_identical(data, sample_names)
end

"""
    pairwise_identical(data::PopData, sample_names::Vector{String})
Return a table of the percent of identical genotypes at each locus
between all pairs of provided individuals.
"""
function pairwise_identical(data::PopData, sample_names::Vector{String})
    errs = ""
    all_samples = samples(data)
    if sample_names != all_samples
        [errs *= "\n  $i" for i in sample_names if i ∉ all_samples]
        errs != "" && error("Samples not found in the PopData: " * errs)
    end
    sample_pairs = pairwise_pairs(sample_names)
    n = length(sample_pairs)
    perc_ident_vec = Vector{Float64}(undef, n)
    n_vec = Vector{Int}(undef, n)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    @inbounds @sync for i in 1:length(sample_pairs)
        Base.Threads.@spawn begin
            @inbounds geno_1 = popdata_idx[(sample_pairs[i][1],)].genotype
            @inbounds geno_2 = popdata_idx[(sample_pairs[i][2],)].genotype
            len_1 = nonmissing(geno_1)
            len_2 = nonmissing(geno_2)
            shared_geno = minimum([len_1, len_2])
                shared_geno = minimum([len_1, len_2])
                shared = sum(skipmissing(geno_1 .== geno_2))
                @inbounds perc_ident_vec[i] = round(shared/shared_geno, digits = 2)
                @inbounds n_vec[i] = shared_geno
            next!(p)
        end
    end
    DataFrame(:sample_1 => map(i -> i[1], sample_pairs), :sample_2 => map(i -> i[2], sample_pairs), :identical => perc_ident_vec, :n => n_vec)
end
