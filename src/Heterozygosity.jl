export ishom, ishet, heterozygosity, het,

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
    return @inbounds all(locus[1] .== locus)
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
    return @inbounds !ishom(locus)
end

@inline function ishet(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishet(locus::T) where T<:GenotypeArray
    return @inbounds ishet.(locus)
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
    1 - sum(allele_freq_vec(data) .^ 2)
end


"""
    het_expected(data::PopData)
    Returns the expected heterozygosity of loci in a PopData object.
"""
function het_expected(data::PopData)
    tmp = @groupby data.loci (:locus, :population) {het_exp = hetero_e(:genotype)}
    @groupby tmp :loci {het_exp = mean(:het_exp |> skipmissing)}
end


"""
    heterozygosity(data::PopData, mode::String = "locus")
Calculate observed and expected heterozygosity in a `PopData` object.
### Modes
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population
## Example
heterozygosity(nancycats(), "population" )
"""
function heterozygosity(data::PopData, mode::String = "locus")
    if mode ∈ ["locus", "loci"]
        tmp = @groupby data.loci (:locus, :population) {het_obs = hetero_o(:genotype), het_exp = hetero_e(:genotype)}
        return @groupby tmp :loci {het_obs = mean(:het_obs |> skipmissing), het_exp = mean(:het_exp |> skipmissing)}

    elseif lowercase(mode) ∈  ["sample", "ind", "individual"]
        return @groupby data.loci :name {het_obs = hetero_o(:genotype)}

    elseif lowercase(mode) ∈  ["pop", "population"]
        #TODO fix this
        return @groupby data.loci :population {het_obs = hetero_o(:genotype), het_exp = hetero_e(allele_freq_vec(:genotype))}
    else
        error("please specify mode \"locus\", \"sample\", or \"population\"")
    end
end

const het = heterozygosity
const He = heterozygosity


"""
    het_sample(data::PopData, individual::String)
Calculate the observed heterozygosity for an individual in a `PopData` object.
Returns an array of heterozygosity values.
"""
function het_sample(data::PopData, individual::String)
    # calculate observed heterozygosity for an individual
    @apply data.loci begin
        @where data.loci :name == individual
        @with {het = hetero_o(:genotype)}
    end
end



const hwe = hwe_test

#TODO

"""
Heterozygosity with uneven population size adjustments
"""
function het_ot(data::PopData)
    tmp = @groupby data.loci (:locus, :population) {het_pop_obs = hetero_o(:genotype)}
    @groupby tmp :locus {het_obs = mean(skipmissing(:het_pop_obs))}
end


function het_ott(data::PopData)
    @groupby data.loci :locus {het_obs = hetero_o(:genotype)}
end
