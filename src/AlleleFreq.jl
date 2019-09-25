"""
  count_alleles(x::PopObj)
Returns an array of `Dicts` of allele counts per locus
"""
function allele_count_beta(x::PopObj)
    y = PopOpt(x)
    tmp = names(y.loci)[1]  # restrict to single locus
    d = Dict()
    for row in y.loci[!, tmp]
        if row === missing
            d[missing] = get!(d, missing, 0) +1
            continue
        else
            for allele in row
                d[allele] = get!(d, allele, 0) +1
            end
        end
    end
    return d
end
