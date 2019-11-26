"""
    geno_freq(locus::Array{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of
genotypes and their frequencies.
"""
function geno_freq(locus::Array{Union{Missing,Tuple},1})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in locus
        if genotype !== missing
        # sum up non-missing genotypes
            d[genotype] = get!(d, genotype, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    geno_freq(locus::SubArray{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in `PopObj` split
by population using `group()`. Returns a `Dict` of genotypes and their frequencies.
"""
function geno_freq(locus::SubArray{Union{Missing,Tuple},1})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in locus
        if genotype !== missing
        # sum up non-missing genotypes
            d[genotype] = get!(d, genotype, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

#=============================================================
REMOVE?

"""
  allele_freq_alpha(locus::PopObj)
Returns an array of `Dicts` of allele counts per locus
"""
function allele_freq_alpha(locus::PopObj)
    y = PopOpt(locus)
    tmp = names(y.loci)[1]  # restrict to single locus for testing
    d = Dict()
    for genotype in y.loci[!, tmp]
        if genotype === missing
            continue
        else
            for allele in genotype
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
    allele_freq_mini(locus::Array{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj`. Returns a `Dict` of
allele's and their frequencies.
"""
function allele_freq_mini(locus::Array{Union{Missing,Tuple},1})
    d = Dict()
    for genotype in locus
        if genotype !== missing
        # sum up alleles
            for allele in genotype
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end


"""
    allele_freq_mini(locus::SubArray{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()`. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq_mini(locus::SubArray{Union{Missing,Tuple},1})
    d = Dict()
    for genotype in locus
        if genotype !== missing
        # sum up alleles
            for allele in genotype
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end
