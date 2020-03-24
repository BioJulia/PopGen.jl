export hwe_test, hwe

"""
    locus_chi_sq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value
"""
function Χ²_locus(locus::T) where T <: AbstractVector
    ## Get expected number of genotypes in a locus
    expected = geno_count_expected(locus)

    ## Get observed number of genotypes in a locus
    observed = geno_count_observed(locus)

    Χ²_stat = expected
    @inbounds for genotype in keys(expected)
        o = get(observed, genotype, 0)
        e = get(expected, genotype, 0)
        Χ²_stat[genotype] = (o - e)^2 / e
    end
    Χ²_stat = values(Χ²_stat) |> sum
    n_geno_exp = length(expected)
    n_alleles = length(unique_alleles(locus))
    #return n_geno_exp, n_alleles
    df = n_geno_exp - n_alleles

    if df > 0
        Χ²_dist = Distributions.Chisq(df)
        p_val = 1 - Distributions.cdf(Χ²_dist, Χ²_stat)
    else
        p_val = missing
    end
    return (Χ²_stat, df, p_val)
end

#TODO
a = @where x.loci :locus == "contig_2784"
a = a.columns.genotype |> collect

"""
    multitest_missing(pvals::Array{Float64,1}, correction::String)
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. Missing values are first removed from the array, the appropriate
correction made, then missing values are re-added to the array at their original
positions. See MultipleTesting.jl docs for full more detailed information.
#### example
`multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")`

### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` or `"b-h"` : Benjamini-Hochberg adjustment
- `"by"` or `"b-y"`: Benjamini-Yekutieli adjustment
- `"bl"` or `"b-l"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forward stop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` or `"b-c"` : Barber-Candès adjustment
"""
@inline function multitest_missing(pvals::Array{<:Union{Missing, <:AbstractFloat},1}, correction::String)
    # get indices of where original missing are
    miss_idx = findall(i -> i === missing, pvals)
    # make seperate array for non-missing P vals
    p_no_miss = skipmissing(pvals) |> collect

    # make a dict of all possible tests and their respective functions
    d = Dict(
        "bonferroni" => Bonferroni(),
        "holm" => Holm(),
        "hochberg" => Hochberg(),
        "bh" => BenjaminiHochberg(),
        "b-h" => BenjaminiHochberg(),
        "by" => BenjaminiYekutieli(),
        "b-y" => BenjaminiYekutieli(),
        "bl" => BenjaminiLiu(),
        "b-l" => BenjaminiLiu(),
        "hommel" => Hommel(),
        "sidak" => Sidak(),
        "forward stop" => ForwardStop(),
        "fs" => ForwardStop(),
        "bc" => BarberCandes(),
        "b-c" => BarberCandes(),
    )

    correct = adjust(p_no_miss, d[lowercase(correction)]) |> Vector{Union{Missing, Float64}}

    # re-add missing to original positions
    @inbounds for i in miss_idx
        @inbounds insert!(correct, i, missing)
    end
    return correct
end


"""
    hwe_test(data::PopObj; by_pop::Bool = false; correction = "none")
Calculate chi-squared test of HWE for each locus and returns observed and
expected heterozygosity with chi-squared, degrees of freedom and p-values for
each locus. Use `by_pop = true` to perform this separately for each population
(default: by_pop = false) and return a NamedTuple of DataFrames with the names
corresponding to the population names. Use `correction =` to specify a P-value
correction method for multiple testing.

#### example
`hwe_test(gulfsharks(), correction = "bh")` \n
`hwe_test(gulfsharks(), by_pop = true, correction = "bh")` \n


### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` or `"b-h"` : Benjamini-Hochberg adjustment
- `"by"` or `"b-y"`: Benjamini-Yekutieli adjustment
- `"bl"` or `"b-l"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forward stop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` or `"b-c"` : Barber-Candès adjustment
"""
@inline function hwe_test(data::PopData; by_pop::Bool = false, correction::String = "none")
    if !by_pop
        tmp = @groupby data.loci :locus {Χ² = Χ²_locus(:genotype)}
        @map tmp {:locus, Χ² = getindex(:Χ², 1), df = getindex(:Χ², 2), P =  getindex(:Χ², 3)}
    else

    end
end

const hwe = hwe_test

## JASON LOOK AT ME PLEEEEEASE

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
