"""
    permute_loci!(data::PopData)
Edits `PopData` in place with loci permuted across populations within
the `.loci` dataframe.
"""
@inline function permute_loci!(data::PopData)
    @inbounds @sync for locus in groupby(data.loci, :locus)
        Base.Threads.@spawn begin
            shuffle!(locus.population)
        end
    end
    data
end

"""
    permute_samples!(data::PopData; meta::Bool = false)
Edits `PopData` in place with samples permuted across populations within
the `.loci` dataframe. Since performance is important for many permutations,
the default is to only edit the `.loci` table in place; use `meta = true`
if you also require the `.meta` dataframe edited in place.
"""
@inline function permute_samples!(data::PopData; meta::Bool = false)
    pops = shuffle(data.meta.population)

    if meta == true
        meta_pops = deepcopy(pops)
        @inbounds for name in groupby(data.meta, :name)
            @inbounds name.population .= pop!(meta_pops)
        end
    end
    @inbounds @sync for name in groupby(data.loci, :name)
        Base.Threads.@spawn begin
            @inbounds name.population .= pop!(pops) 
        end
    end
    data
end


@inline function permute_samples!(data::AbstractDataFrame, popnames::Vector{String})
    pops = shuffle(popnames)

    @inbounds @sync for name in groupby(data, :name)
        Base.Threads.@spawn begin 
            @inbounds name.population .= pop!(pops)
        end
    end
    data
end


"""
    permute_genotypes!(data::PopData; by::String = "locus")
Edits `PopData` in place with genotypes permuted across individuals within
the `.loci` dataframe. Use `by = "population"` (or `"pop"`) to permute genotypes
within populations.
"""
@inline function permute_genotypes!(data::PopData; by::String = "locus")
    #establish mode of operation
    if by in ["locus", "loci"]
        groupings = :locus
    elseif by in ["population", "pop"]
        groupings = [:locus, :population]
    else
        error("Please choose from either \"locus\" or \"population\" run methods.")
    end
    @inbounds @sync for grp in groupby(data.loci, groupings)
        Base.Threads.@spawn begin
            grp.genotype .= strict_shuffle!(grp.genotype)
        end
    end
    data
end


"""
    permute_alleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
Edits `PopData` in place with alleles permuted and reconstructed into genotypes
for each locus within the `.loci` dataframe. Use `by = "population"` (or `"pop"`)
to permute alleles within populations. If `ploidy` is not provided (default `ploidy = nothing`),
then ploidy will be identified from the PopData. If performance is important,
it would be best to identify ploidy in advance and set it to a specific integer.
"""
@inline function permute_alleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
    #establish mode of operation
    if by in ["locus", "loci"]
        groupings = :locus
    elseif by in ["population", "pop"]
        groupings = [:locus, :population]
    else
        error("Please choose from either \"locus\" or \"population\" run methods.")
    end

    if ploidy == nothing
        tmp = unique(data.meta.ploidy)
        length(tmp) > 1 && error("This permutation method is not appropriate for mixed-ploidy data")
        ploidy = tmp[1]
    end

    @inbounds @sync for grp in groupby(data.loci, groupings)
        Base.Threads.@spawn begin
            alle = shuffle(alleles(grp.genotype))
            new_genos = Tuple.(Base.Iterators.partition(alle, ploidy))
            (@view grp.genotype[.!ismissing.(grp.genotype)]) .= new_genos
        end
    end
    data
end
