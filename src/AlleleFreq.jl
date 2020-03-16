"""
    alleles(locus::T) where T<:AbstractArray
Return an array of all the non-missing alleles of a locus.
"""
@inline function alleles(locus::T) where T<:AbstractArray
    Base.Iterators.flatten(locus |> skipmissing) |> collect
end


"""
    allele_freq(locus::T) where T<:AbstractArray
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object.
"""
@inline function allele_freq(locus::T) where T<:AbstractArray
    d = Dict{Int16,Float32}()
    flat_alleles = alleles(locus)
    len = length(flat_alleles)
    @inbounds Base.@async for allele in flat_alleles
        @inbounds d[allele] = get!(d, allele, 0) + 1/len
    end
    return d
end


@inline function allele_freq_vec(locus::T) where T<:AbstractArray
        flat_alleles = alleles(locus)
        len = length(flat_alleles)
        d = [count(i -> i == j, flat_alleles) for j in unique(flat_alleles)]
        return d ./ len
end

"""
    allele_freq(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of allele frequencies
of that locus per population.
### Example
```
cats = nancycats()
allele_freq(cats, "fca8")
allele_freq(cats, "fca8", population = true)
```
"""
function allele_freq(data::PopData, locus::String, population::Bool=false)
    if !population
        @apply data.loci begin
            @where :locus == locus
            @with allele_freq(:genotype)
        end
    else
        @apply data.loci :population begin
            @where :locus == locus
            @with {freq = allele_freq(:genotype)}
        end
    end
end


"""
    geno_freq(locus::T) where T<:AbstractArray
Return a `Dict` of genotype frequencies of a single locus in a
`PopData` object.
"""
@inline function geno_freq(locus::T) where T<:AbstractArray
    # conditional testing if all genos are missing
    all(ismissing.(locus)) && return missing

    d = Dict{Tuple, Float32}()
    @inbounds for genotype in skipmissing(locus)
        # sum up non-missing genotypes
        d[genotype] = get!(d, genotype, 0) + 1
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end


"""
    geno_freq(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.
### Example
```
cats = nancycats()
geno_freq(cats, "fca8")
geno_freq(cats, "fca8", population = true)
```
"""
function geno_freq(data::PopData, locus::String, population::Bool=false)
    if !population
        @apply data.loci begin
            @where :locus == locus
            @with geno_freq(:genotype)
        end
    else
        @apply data.loci :population begin
            @where :locus == locus
            @with {freq = geno_freq(:genotype)}
        end
    end
end


"""
    get_genotype(data::PopObj; sample::String, locus::String)
Return the genotype of a specific sample for specific locus in a `PopData` object.

`get_genotype(nancycats, sample = "N115" , locus = "fca8")`
"""
function get_genotype(data::PopData; sample::String, locus::String)
    @where data.loci :name == sample && :locus == locus
end

"""
    get_sample_genotypes(data::PopData, sample::String)
Return all the genotypes of a specific sample in a `PopObj`.
This is an extension for the internal function `get_genotype`.
```
cats = nancycats()
get_sample_genotypes(cats, "N115")
```
"""
function get_sample_genotypes(data::PopObj, sample::String)
    @where data.loci :name == sample
end
