"""
  count_alleles(x::PopObj)
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
end
