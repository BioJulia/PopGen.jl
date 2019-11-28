"""
    geno_freq(locus::Array{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of
genotypes and their frequencies.
"""
function geno_freq(locus::Array{Union{Missing,Tuple},1})
    d = Dict()
    # conditional testing if all genos are missing
    all(ismissing.(locus)) == true && return missing
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
    all(ismissing.(locus)) == true && return missing
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
    allele_freq_mini(locus::Array{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj`. Returns a `Dict` of
allele's and their frequencies.
"""
function allele_freq_mini(locus::Array{Union{Missing,Tuple},1})
    d = Dict()
    for genotype in skipmissing(locus)
    # sum up alleles
        for allele in genotype
            d[allele] = get!(d, allele, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end

#=================================
function allele_freq_micro(locus::Array{Union{Missing,Tuple},1})
    nonmissing = skipmissing(locus)
    flattened_alleles = collect(Iterators.flatten(nonmissing))
    # sum up alleles
    count_d = countmap(flattened_alleles)
    total = values(count_d) |> sum    # sum of all non-missing alleles
    out_d = Dict()    # the dict of freqs that will be returned
    [out_d[i] = count_d[i] / total for i in keys(count_d)]  # allele count / sum
    return out_d
end
===================================#

"""
    allele_freq_micro(locus::SubArray{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()`. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq_mini(locus::SubArray{Union{Missing,Tuple},1})
    d = Dict()
    for genotype in skipmissing(locus)
    # sum up alleles
        for allele in genotype
            d[allele] = get!(d, allele, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end
