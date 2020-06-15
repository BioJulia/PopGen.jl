
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
    permute_genotypes!(data::PopData)
Edits `PopData` in place with genotypes permuted across individuals within
the `.loci` dataframe.
"""
@inline function permute_genotypes!(data::PopData)
    @inbounds for locus in groupby(data.loci, :locus)
        locus.genotype .= strict_shuffle!(locus.genotype)
    end
end


"""
    permute_genotypes_pop!(data::PopData)
Edits `PopData` in place with genotypes permuted across individuals within
populations within the `.loci` dataframe.
"""
@inline function permute_genotypes_pop!(data::PopData)
    @inbounds for locus in groupby(data.loci, [:locus, :population])
        locus.genotype .= strict_shuffle!(locus.genotype)
    end
end

#WIP
"""
    permute_alleles!(data::PopData)

"""
@inline function permute_alleles!(data::PopData)
    @inbounds for locus in groupby(data.loci, :locus)
        alle = alleles(locus.genotype, miss = true)
        return shuffle(alle)
    end
end
