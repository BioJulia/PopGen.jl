"""
    geno_freq(x::Array{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in a `PopObj`
"""
function geno_freq(x::Array{Union{Missing,Tuple},1})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(x) |> unique == [true] && return missing
    for row in x
        # sum up missing
        if row === missing
            continue
        else
        # sum up non-missing genotypes
            d[row] = get!(d, row, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    geno_freq(x::SubArray{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in `PopObj` split
by population using group()
"""
function geno_freq(x::SubArray{Union{Missing,Tuple},1})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(x) |> unique == [true] && return missing
    for row in x
        # sum up missing
        if row === missing
            continue
        else
        # sum up non-missing genotypes
            d[row] = get!(d, row, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

#=============================================================
REMOVE?

"""
  allele_freq_alpha(x::PopObj)
Returns an array of `Dicts` of allele counts per locus
"""
function allele_freq_alpha(x::PopObj)
    y = PopOpt(x)
    tmp = names(y.loci)[1]  # restrict to single locus for testing
    d = Dict()
    for row in y.loci[!, tmp]
        if row === missing
            continue
        else
            for allele in row
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum
    [d[i] = d[i] / total for i in keys(d)]
    return d
end
==============================================================#

"""
    allele_freq_mini(x::Array{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj`
"""
function allele_freq_mini(x::Array{Union{Missing,Tuple},1})
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
