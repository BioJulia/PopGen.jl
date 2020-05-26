export ishom, ishet, heterozygosity, het

"""
```
ishom(locus::T) where T <: GenotypeArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't, and `missing` if it's `missing`. The vector version
simply broadcasts the function over the elements.
"""
@inline function ishom(locus::Genotype)
    # if allele 1 equels all others, return true
    return all(@inbounds locus[1] .== locus)
end

@inline function ishom(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishom(locus::T) where T<:GenotypeArray
    return @inbounds ishom.(locus)
end


"""
```
ishet(locus::T) where T <: GenotypeArray
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

@inline function ishet(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishet(locus::T) where T<:GenotypeArray
    return @inbounds ishet.(locus)
end

#TODO add to docs (API)

"""
    gene_diversity_nei87(het_exp::Union{Missing,AbstractFloat}, het_obs::Union{Missing,AbstractFloat}, n::U where U<:Signed)
Calculate overall gene diversity with the adjustment/correction
given by Nei:

Hₜ = 1 −sum(pbar²ᵢ + Hₛ/(ñ * np) − Het_obs/(2ñ*np))
- _ñ_ is the number of genotypes for a locus for a population
- _np_ is the number of genotypes of a locus across all populations
    - i.e. sum(_ñ_)
- _pbar²_ is the observed homozygosity of a locus for that population
- _Hₛ_ is the within population gene diversity given by:
    - Hₛ =  ñ/(ñ-1) * (1 - sum(pbar²ᵢ - Het_observed / 2ñ))

(Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press).
"""
function gene_diversity_nei87(het_exp::Union{Missing,AbstractFloat}, het_obs::Union{Missing,AbstractFloat}, n::U where U<:Signed)
    het_exp - (het_obs/n/2) * n/(n-1)
end

"""
    hetero_o(data::T) where T <: GenotypeArray
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined
as genotypes returning `true` for `ishet()`. This is numerically feasible because
`true` values are mathematically represented as `1`, whereas `false` are represented
as `0`.
"""
@inline function hetero_o(data::T) where T <: GenotypeArray
    adjusted_vector = ishet(data) |> skipmissing
    isempty(adjusted_vector) ? missing : mean(adjusted_vector)
end


"""
    hetero_e(allele_freqs::Vector{T}) where T <: GenotypeArray
Returns the expected heterozygosity of an array of genotypes,
calculated as 1 - sum of the squared allele frequencies.
"""
function hetero_e(data::T) where T <: GenotypeArray
    1 - sum(@inbounds @avx allele_freq_vec(data) .^ 2)
end


"""
    heterozygosity(data::PopData, mode::String = "locus")
Calculate observed and expected heterozygosity in a `PopData` object. For loci,
heterozygosity is calculated in the Nei fashion, such that heterozygosity is
calculated as the average over heterozygosity per locus per population.
### Modes
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population
## Example
heterozygosity(nancycats(), "population" )
"""
function heterozygosity(data::PopData, mode::String = "locus")
    if mode ∈ ["locus", "loci"]
        tmp = @groupby data.loci (:locus, :population) {n_tmp = nonmissing(:genotype), het_pop_obs = hetero_o(:genotype), het_pop_exp = hetero_e(:genotype)}
        @groupby tmp :locus {n = sum(:n_tmp), het_obs = mean(skipmissing(:het_pop_obs)), het_exp = mean(skipmissing(:het_pop_exp))}

    elseif lowercase(mode) ∈  ["sample", "ind", "individual"]
        return @groupby data.loci :name {n = nonmissing(:genotype), het_obs = hetero_o(:genotype)}

    elseif lowercase(mode) ∈  ["pop", "population"]
        return @groupby data.loci :population {n = nonmissing(:genotype), het_obs = hetero_o(:genotype), het_exp = hetero_e(:genotype)}

    else
        error("please specify mode \"locus\", \"sample\", or \"population\"")
    end
end

const het = heterozygosity

"""
    het_sample(data::PopData, individual::String)
Calculate the observed heterozygosity for an individual in a `PopData` object.
Returns an array of heterozygosity values.
"""
function het_sample(data::PopData, individual::String)
    data.loci[data.loci.name .== individual, :genotype] |> hetero_o
end
