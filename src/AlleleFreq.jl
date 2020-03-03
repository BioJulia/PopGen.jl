"""
    alleles(locus::T) where T<:AbstractArray
Return an array of all the non-missing alleles of a locus.
"""
@inline function alleles(locus::T) where T<:AbstractArray
    Base.Iterators.flatten(locus |> skipmissing) |> collect
end


"""
    allele_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate allele counts for a single locus of a `PopObj`. Returns a
 `Dict` of allele's and their frequencies.
"""
@inline function allele_freq(locus::T) where T<:AbstractArray
    d = Dict{String,Float32}()
    flat_alleles = alleles(locus)
    uniq_alleles = unique(flat_alleles)
    allele_counts = [count(i -> i == j, flat_alleles) for j in uniq_alleles]
    freq = allele_counts ./ sum(allele_counts)
    for (i,j) in enumerate(uniq_alleles)
        d[string(j)] = freq[i]
    end
    return d
end

#TODO
"""
    allele_matrix(data::PopObj, allele::Int)
Convert a DataFrame of genotypes from a PopObj into a matrix of a single
allele. Use `allele` to specify which allele you want isolated. Intended
to be used to create separate allele matrices for comparing alleles 1 and
2, such as assessing heterozygosity.
"""
function allele_matrix(data::T, allele::Int) where T <: DataFrame
    geno_matrix = Matrix(data)
    map(geno_matrix) do i
        i === missing && return missing
        i[allele]
    end
end


"""
    geno_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of
genotypes and their frequencies.
"""
@inline function geno_freq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    # conditional testing if all genos are missing
    all(ismissing.(locus)) == true && return missing
    d = Dict()
    @inbounds for genotype in skipmissing(locus)
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
@inline function geno_freq(locus::Vector{<:Union{Missing, Tuple{Vararg}}})
    # conditional testing if all genos are missing
    all(ismissing.(locus)) == true && return missing
    d = Dict()
    @inbounds for genotype in skipmissing(locus)
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
@inline function geno_freq(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    # conditional testing if all genos are missing
    all(ismissing.(locus)) == true && return missing
    d = Dict()
    @inbounds for genotype in skipmissing(locus)
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
@inline function geno_freq(locus::SubArray{<:Union{Missing, Tuple{Vararg}}})
    # conditional testing if all genos are missing
    all(ismissing.(locus)) == true && return missing
    d = Dict()
    @inbounds for genotype in skipmissing(locus)
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
