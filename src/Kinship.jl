"""
    allele_freq(x::Tuple)
Calculate allele frequency for a single locus of a single sample. Returns a
`Dict` of alleles and their frequencies.
"""
function allele_freq(x::Array{Union{Missing,Tuple},1})
    d = Dict()
    for row in x
        # sum up missing
        if row === missing
            continue
        else
        # sum up alleles
            for allele in row
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end
