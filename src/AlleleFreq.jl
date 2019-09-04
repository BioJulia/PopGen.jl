"""
  count_alleles(x::PopObj)
Returns an array of `Dicts` of genotype counts per locus
"""
function count_alleles(x::PopObj)
    tmp = []
    for locus in x.loci
        y = x.genotypes[locus]
        d=Dict([("$i",count(x->x==i,y)) for i in unique(y)])
        push!(tmp, d)
    end
    return tmp
end
