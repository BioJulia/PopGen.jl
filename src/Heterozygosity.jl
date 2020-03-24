export ishom, ishet, heterozygosity, het,

"""
```
ishom(locus::T) where T<:AbstractVector
ishom(locus::NTuple{N,T} where N where T <: Signed)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't. The vector version simply broadcasts the function over the
elements.
"""
@inline function ishom(locus::NTuple{N,T} where N where T <: Signed)
    # if allele 1 equels all others, return true
    return @inbounds all(locus[1] .== locus)
end

@inline function ishom(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishom(locus::T) where T<:AbstractVector
    return @inbounds ishom.(locus)
end


"""
```
ishet(locus::T) where T<:AbstractVector
ishet(locus::NTuple{N,T} where N where T <: Signed)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if
it is, `false` if it isn't. The vector version simply broadcasts the function over the
elements. Under the hood, this function is simply `!ishom`.
"""
@inline function ishet(locus::NTuple{N,T} where N where T <: Signed)
    return @inbounds !ishom(locus)
end

@inline function ishet(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishet(locus::T) where T<:AbstractVector
    return @inbounds ishet.(locus)
end


"""
    hetero_o(data::T) where T <: AbstractVector
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined
as genotypes returning `true` for `ishet()`. This is numerically feasible because
`true` values are mathematically represented as `1`, whereas `false` are represented
as `0`.
"""
@inline function hetero_o(data::T) where T <: AbstractVector
    adjusted_vector = ishet(data) |> skipmissing
    isempty(adjusted_vector) ? missing : mean(adjusted_vector)
end


"""
    hetero_e(allele_freqs::Vector{T}) where T <: AbstractFloat
Returns the expected heterozygosity of an array of genotypes,
calculated as 1 - sum of the squared allele frequencies.
"""
function hetero_e(data::T) where T <: AbstractVector
    1 - sum(allele_freq_vec(data) .^ 2)
end


"""
    het_expected(data::PopData)
    Returns the expected heterozygosity of loci in a PopData object.
"""
function het_expected(data::PopData)
    @groupby data.loci :locus {het_exp = hetero_e(:genotype)}
end


"""
    heterozygosity(data::PopData, mode::String = "locus")
Calculate observed and expected heterozygosity in a `PopData` object.
## Example
heterozygosity(nancycats(), "population" )
### Modes
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population
"""
function heterozygosity(data::PopData, mode::String = "locus")
    if mode ∈ ["locus", "loci"]
        return @groupby data.loci :locus {het_obs = hetero_o(:genotype), het_exp = hetero_e(:genotype)}

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
Calculate the observed heterozygosity for an individual in a `PopObj`. Returns
an array of heterozygosity values.
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

function het_obstest(data::PopObj)
    het_vals = Vector{Float64}()
    for locus in eachcol(data.loci, false)
        # get genotype freqs at locus
        a = geno_freq(locus)
        tmp = 0
        # total number of non-missing samples per locus
        n = skipmissing(locus) |> collect |> length
        genos = keys(a) |> collect
        for geno in genos
            if length(unique(geno)) != 1
                # if true, add freq to total in tmp
                tmp += a[geno]
            end
        end
        # include het correction
        het_adjust = tmp * (n/ (n-1))
        push!(het_vals, het_adjust)
    end

    return het_vals
end


@inline function hwe_test(data::PopData; by_pop::Bool = false, correction::String = "none")
    if !by_pop
        tmp = @groupby data.loci :locus {Χ² = locus_chi_sq(:genotype)}
        @map tmp {:locus, Χ² = getindex(:Χ²,1), df = getindex(:Χ²,2), P =  getindex(:Χ²,3)}
    else

    end
end

jdb = hwe_test(x)
df = by(y, :locus, :genotype => locus_chi_sq)
