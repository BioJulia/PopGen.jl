
"""
    permute_loci!(data::PopData)
Edits `PopData` in place with loci permuted across populations within
the `.loci` dataframe.
"""
@inline function permute_loci!(data::PopData)
    @inbounds for locus in groupby(data.loci, :locus)
        shuffle!(locus.population)
    end
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

    @inbounds for name in groupby(data.loci, :name)
        @inbounds name.population .= pop!(pops)
    end

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
    else
        groupings = [:locus, :population]
    end
    @inbounds for locus in groupby(data.loci, groupings)
        locus.genotype .= strict_shuffle!(locus.genotype)
    end
end


#WIP
"""
    permute_alleles!(data::PopData)

"""
@inline function permute_alleles!(data::PopData; ploidy::Int = 2, by::String = "locus")
    #establish mode of operation
    if by in ["locus", "loci"]
        groupings = :locus
    else
        groupings = [:locus, :population]
    end

    @inbounds for locus in groupby(data.loci, groupings)
        alle = shuffle(alleles(locus.genotype))
        new_genos = Tuple.(Base.Iterators.partition(alle, ploidy))
        @inbounds locus.genotype[@view .!ismissing.(locus.genotype)] .= new_genos
    end
end
