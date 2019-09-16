"""
  count_alleles(x::PopObj)
Returns an array of `Dicts` of allele counts per locus
"""
function count_alleles(x::PopObj)
    al_counts = []
    for loci in collect(keys(x.genotypes))
        d = Dict()
        for locus in x.genotypes[loci]
            for allele in locus
                d[allele] = get!(d, allele, 0) + 1
            end
        end
        push!(al_counts, d)
    end
    return al_counts
end
