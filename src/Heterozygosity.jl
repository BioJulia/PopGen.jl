export ishom, ishet, heterozygosity, het


#TODO how to treat haploids?
"""
```
ishom(locus::T) where T <: GenoArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't, and `missing` if it's `missing`. The vector version
simply maps the function over the elements.
"""
@inline function ishom(locus::Genotype)
    # if the first equals all others, return true
    return all(@inbounds first(locus) .== locus)
end

ishom(locus::Missing) = missing

@inline function ishom(locus::T) where T<:GenoArray
    return @inbounds map(ishom, locus)
end

@inline function ishom(locus::T) where T<:Base.SkipMissing
    return @inbounds map(ishom, locus)
end


"""
    ishom(locus::Genotype, allele::Signed)
    ishom(loci::GenoArray, allele::Signed)
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele`.
"""
function ishom(geno::T, allele::U) where T<:Genotype where U<:Integer
    ∈(allele, geno) & ishom(geno) ? true : false
end

ishom(geno::T, allele::U) where T<:GenoArray where U<:Integer = map(i -> ishom(i, allele), geno)

ishom(geno::Missing, allele::U) where U<:Integer = missing

"""
```
ishet(locus::T) where T <: GenoArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if
it is, `false` if it isn't. The vector version simply broadcasts the function over the
elements. Under the hood, this function is simply `!ishom`.
"""
@inline function ishet(locus::Genotype)
    return !ishom(locus)
end

ishet(locus::Missing) = missing


@inline function ishet(locus::T) where T<:GenoArray
    return @inbounds map(ishet, locus)
end


@inline function ishet(locus::T) where T<:Base.SkipMissing
    return @inbounds map(ishet, locus)
end


"""
    ishet(locus::Genotype, allele::Signed)
    ishet(loci::GenoArray, allele::Signed)
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele`. 
"""
function ishet(geno::T, allele::U) where T<:Genotype where U<:Integer
    ∈(allele, geno) & !ishom(geno) ? true : false
end

ishet(geno::T, allele::U) where T<:GenoArray where U<:Integer = map(i -> ishet(i, allele), geno)

ishet(geno::Missing, allele::U) where U<:Integer = missing


"""
    counthet(geno::T, allele::Int) where T<:GenoArray
    counthet(geno::T, allele::AbstractVector{U}) where T<:GenoArray where U<:Integer
Given a `GenoArray`, count the number of times `allele` appears in the
heterozygous state.
"""
function counthet(geno::T, allele::U) where T<:GenoArray where U<:Integer
    mapreduce(i -> ishet(i, allele), +, skipmissing(geno))
end


function counthet(geno::T, allele::AbstractVector{U}) where T<:GenoArray where U<:Integer
    tmp = skipmissing(geno)
    isempty(tmp) && return fill(0, length(allele))
    map(_allele -> mapreduce(i -> ishet(i, _allele), +, tmp), allele)
end


"""
    counthom(geno::T, allele::Int) where T<:GenoArray
Given a `GenoArray`, count the number of times `allele` appears in the
homozygous state.
"""
function counthom(geno::T, allele::U) where T<:GenoArray where U <: Signed
    mapreduce(i -> ishom(i, allele), +, skipmissing(geno))
end

function counthom(geno::T, allele::AbstractVector{U}) where T<:GenoArray where U<:Integer
    tmp = skipmissing(geno)
    isempty(tmp) && return fill(0, length(allele))
    map(_allele -> mapreduce(i -> ishom(i, _allele), +, tmp), allele)
end

"""
    gene_diversity_nei87(het_exp::Union{Missing,AbstractFloat}, het_obs::Union{Missing,AbstractFloat}, n::Union{Integer, Float64}, corr::Bool = true)
Calculate overall gene diversity with the adjustment/correction given by Nei:

Hₜ = 1 −sum(pbar²ᵢ + Hₛ/(ñ * np) − Het_obs/(2ñ*np))
- _ñ_ is the number of genotypes for a locus for a population
- _np_ is the number of genotypes of a locus across all populations
    - i.e. sum(_ñ_)
- _pbar²_ is the observed homozygosity of a locus for that population
- _Hₛ_ is the within population gene diversity given by:
    - Hₛ =  ñ/(ñ-1) * (1 - sum(pbar²ᵢ - Het_observed / 2ñ))

(Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press).
use `corr = false` to ignore sample-size correction `* n/(n-1)`.
"""
@inline function gene_diversity_nei87(het_exp::T, het_obs::T, n::Union{Integer,T}; corr::Bool = true) where T<: AbstractFloat
    corr_val = corr ? n/(n-1.0) : 1.0
    return @fastmath (het_exp - (het_obs/n/2.0)) * corr_val
end

@inline function gene_diversity_nei87(het_exp::AbstractFloat, het_obs::Missing, n::Union{Integer,AbstractFloat}; corr::Bool = true)
    return missing
end

@inline function gene_diversity_nei87(het_exp::Missing, het_obs::AbstractFloat, n::Union{Integer,AbstractFloat}; corr::Bool = true)
    return missing
end

@inline function gene_diversity_nei87(het_exp::Missing, het_obs::Missing, n::Union{Integer,AbstractFloat}; corr::Bool = true)
    return missing
end

"""
    hetero_o(data::T) where T <: GenoArray
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined
as genotypes returning `true` for `ishet()`. This is numerically feasible because
`true` values are mathematically represented as `1`, whereas `false` are represented
as `0`.
"""
@inline function hetero_o(data::T) where T <: GenoArray
    adjusted_vector = ishet(data) |> skipmissing
    isempty(adjusted_vector) ? missing : mean(adjusted_vector)
end


"""
    hetero_e(allele_freqs::Vector{T}) where T <: GenoArray
Returns the expected heterozygosity of an array of genotypes,
calculated as 1 - sum of the squared allele frequencies.
"""
@inline function hetero_e(data::T) where T <: GenoArray
    all(ismissing.(data)) == true ? missing : 1.0 - mapreduce(i -> i^2, + , allele_freq_vec(data))
end


"""
    heterozygosity(data::PopData, by::String = "locus")
Calculate observed and expected heterozygosity in a `PopData` object. For loci,
heterozygosity is calculated in the Nei fashion, such that heterozygosity is
calculated as the average over heterozygosity per locus per population.
### Modes
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population
## Example
heterozygosity(@nancycats, "population" )
"""
function heterozygosity(data::PopData, by::String = "locus")
    if by ∈ ["locus", "loci"]
        tmp = DataFrames.combine(
                groupby(data.loci, [:locus, :population]),
                :genotype => nonmissing => :n_tmp,
                :genotype => hetero_o => :het_pop_obs,
                :genotype => hetero_e => :het_pop_exp
            )
        return DataFrames.combine(
                groupby(tmp, :locus),
                :n_tmp => sum => :n,
                :het_pop_obs => (h_o -> mean(skipmissing(h_o))) => :het_obs,
                :het_pop_exp => (h_e -> mean(skipmissing(h_e))) => :het_exp
            )

    elseif lowercase(by) ∈  ["sample", "ind", "individual"]
        return DataFrames.combine(
                groupby(data.loci, :name),
                :genotype => nonmissing => :n,
                :genotype => hetero_o => :het_obs
            )

    elseif lowercase(by) ∈  ["pop", "population"]
        return DataFrames.combine(
                groupby(data.loci, :population),
                :genotype => nonmissing => :n,
                :genotype => hetero_o => :het_obs,
                :genotype => hetero_e => :het_exp
            )
    else
        error("please specify by = \"locus\", \"sample\", or \"population\"")
    end
end

const het = heterozygosity


#NOTE this is not intended to be performant. It's a convenience function. 
"""
    het_sample(data::PopData, individual::String)
Calculate the observed heterozygosity for an individual in a `PopData` object.
"""
@inline function het_sample(data::PopData, individual::String)
    data.loci[data.loci.name .== individual, :genotype] |> hetero_o
end
