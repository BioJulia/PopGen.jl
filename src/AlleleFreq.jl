"""
    allele_avg(data::PopObj, rounding::Bool = true)
Returns a NamedTuple of the average number of alleles ('avg') and standard
deviation (`stdev`) of a `PopObj`. Use `false` as second argument (no keyword)
to not round results. Default (`true`) rounds to 4 digits.
"""
function allele_avg(data::PopObj, rounding::Bool = true)
    all_dicts = map(allele_freq, eachcol(data.loci))
    just_alleles = map(i -> keys(i) |> collect, all_dicts)
    num_alleles = map(length, just_alleles)
    avg = mean(num_alleles)
    sd = std(num_alleles)
    if rounding == true
        return (avg = round(avg, digits = 4), stdev = round(sd, digits = 4))
    else
        return (avg = avg, stdev = sd)
    end
end

"""
    allele_freq(genotype::Union{Missing,Tuple})
Calculate allele frequency for a single locus of a single sample. Returns a
`Dict` of alleles and their frequencies.
"""
function allele_freq(genotype::Union{Missing,Tuple})
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
    allele_freq(locus::Vector{Union{Missing, Tuple}})
Calculate allele counts for a single locus of a `PopObj`. Returns a `Dict` of
allele's and their frequencies.
"""
function allele_freq(locus::Array{Union{Missing,Tuple},1})
    d = Dict()
    for geno in locus
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
    allele_freq(locus::SubArray{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj` split by population
using `group()`. Returns a `Dict` of allele's and their frequencies.
"""
function allele_freq(locus::SubArray{Union{Missing,Tuple},1})
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
    geno_freq(locus::Vector{Union{Missing, Tuple}})
Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of
genotypes and their frequencies.
"""
function geno_freq(locus::Vector{Union{Missing,Tuple}})
    d = Dict()
    # conditional testing if all genos are missing
    ismissing.(locus) |> unique == [true] && return missing
    for genotype in locus
        genotype === missing && continue

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
        genotype === missing && continue

        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
        end
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
