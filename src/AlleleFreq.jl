"""
  count_alleles(x::PopObj)
<<<<<<< HEAD
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
=======
Returns an array of `Dicts` of genotype counts per locus
"""
function count_genotypes(x::PopObj)
    tmp = []
    d = Dict()
    for locus in x.loci
        y = x.genotypes[locus]
        d[locus] = [(genotype = i, count=count(x->x==i,y)) for i in unique(y)]
    end
    return d
end

function count_genotypes_df(x::PopObj)
    loci = []
    genotypes = []
    counts = []
    for locus in x.loci
        y = x.genotypes[locus]
        append!(genotypes, unique(y))
        append!(loci, fill(locus, length(unique(y))))
        append!(counts, [count(x->x==i,y) for i in unique(y)])
    end
    return DataFrame(locus = string.(loci), genotype = genotypes, count = Int.(counts))
>>>>>>> master
end
