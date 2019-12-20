"""
    allele_freq(genotype::Union{Missing, NTuple{N,<:Integer}}) where N
Calculate allele frequency for a single locus of a single sample. Returns a
`Dict` of alleles and their frequencies.
"""
function allele_freq(genotype::Union{Missing, NTuple{N,<:Integer}}) where N
    x === missing && return missing
    d = Dict()
    for allele in genotype
        # sum up alleles
        d[allele] = get!(d, allele, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end

"""
    allele_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate allele counts for a single locus of a `PopObj`. Returns a `Dict` of
allele's and their frequencies.
"""
function allele_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
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

"""
    allele_freq(locus::Vector{<:Union{Missing, Tuple{Vararg}}}
Calculate allele counts for a single locus of a `PopObj` with unequal ploidy across
samples. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq(locus::Vector{<:Union{Missing, Tuple{Vararg}}})
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

"""
    allele_freq(locus::SubArray{<:Union{Missing,NTuple{N,<:Integer}}}) where N
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()`. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq(locus::SubArray{<:Union{Missing,NTuple{N,<:Integer}}}) where N
    d = Dict()
    for genotype in locus
        genotype === missing && continue
        # sum up alleles
        for allele in genotype
            d[allele] = get!(d, allele, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end

"""
    allele_freq(locus::SubArray{<:Union{Missing,Tuple{Vararg}}})
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()` with unequal ploidy across samples. Returns a `Dict` of
allele's and their frequencies.
"""
function allele_freq(locus::SubArray{<:Union{Missing,Tuple{Vararg}}})
    d = Dict()
    for genotype in locus
        genotype === missing && continue
        # sum up alleles
        for allele in genotype
            d[allele] = get!(d, allele, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end

"""
    geno_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of
genotypes and their frequencies.
"""
function geno_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in skipmissing(locus)
        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    geno_freq(locus::Vector{<:Union{Missing, Tuple{Vararg}}})
Calculate genotype frequencies of all loci in a `PopObj` with unequal ploidy across
samples. Returns a `Dict` of genotypes and their frequencies.
"""
function geno_freq(locus::Vector{<:Union{Missing, Tuple{Vararg}}})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in skipmissing(locus)
        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    geno_freq(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate genotype frequencies of all loci in `PopObj` split
by population using `group()`. Returns a `Dict` of genotypes and their frequencies.
"""
function geno_freq(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in skipmissing(locus)
        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    geno_freq(locus::SubArray{<:Union{Missing, Tuple{Vararg}}})
Calculate genotype frequencies of all loci in `PopObj` split by population
using `group()` with unequal ploidy across samples. Returns a `Dict` of
genotypes and their frequencies.
"""
function geno_freq(locus::SubArray{<:Union{Missing, Tuple{Vararg}}})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in skipmissing(locus)
        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end

"""
    get_genotype(data::PopObj; sample::String, locus::String)
Return the genotype of a specific sample for specific locus in a `PopObj`.
This is a barebones variant to `isolate_genotypes` that is ~1000x faster and
suitable for development.

`get_genotype(nancycats, sample = "N115" , locus = "fca8")`

"""
function get_genotype(data::PopObj; sample::String, locus::String)
    idx = findfirst(i -> i == sample, data.samples.name)
    return getindex(data.loci[!, Symbol(locus)], idx)
end

"""
    get_sample_genotypes(data::PopObj, sample::String)
Return all the genotypes of a specific sample in a `PopObj`.
This is an extension for the internal function `get_genotype`.

`get_sample_genotypes(nancycats(), "N115")`
"""
function get_sample_genotypes(data::PopObj, sample::String)
    idx = findfirst(i -> i == sample, data.samples.name)
    map(i -> i[idx] , eachcol(data.loci))
end
