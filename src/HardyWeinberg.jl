export hwe_test, hwe

"""
    locus_chi_sq(locus::T) where T <: GenoArray
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value
"""
function chisq_locus(locus::T) where T <: GenoArray
    ## Get expected number of genotypes in a locus
    expected = geno_count_expected(locus)

    ## Get observed number of genotypes in a locus
    observed = geno_count_observed(locus)

    chisq_stat = expected
    @inbounds for genotype in keys(expected)
        o = get(observed, genotype, 0)
        e = get(expected, genotype, 0)
        chisq_stat[genotype] = (o - e)^2 / e
    end
    chisq_stat = values(chisq_stat) |> sum
    n_geno_exp = length(expected)
    n_alleles = length(unique_alleles(locus))
    df = n_geno_exp - n_alleles

    if df > 0
        chisq_dist = Distributions.Chisq(df)
        p_val = 1.0 - Distributions.cdf(chisq_dist, chisq_stat)
    else
        p_val = missing
    end
    return (chisq_stat, df, p_val)
end


"""
    hwe_test(data::PopData; by::String = "locus"; correction = "none")
Calculate chi-squared test of HWE for each locus and returns observed and
expected heterozygosity with chi-squared, degrees of freedom and p-values for
each locus. Use `by = "population"` to perform this separately for each population
(default: `by = "locus"`). Use `correction =` to specify a P-value
adjustment method for multiple testing.

#### example
`hwe_test(gulfsharks(), correction = "bh")` \n
`hwe_test(gulfsharks(), by = "population", correction = "bh")` \n


### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"`  : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
"""
@inline function hwe_test(data::PopData; by::String = "locus", correction::String = "none")
    if by == "locus"
        tmp =DataFrames.combine(
            groupby(data.loci, :locus),
            :genotype => chisq_locus => :chisq
        )
        out_table = select(tmp, :locus, 
            :chisq => (i -> getindex.(i, 1)) => :chisq, 
            :chisq => (i -> getindex.(i, 2)) => :df, 
            :chisq => (i -> getindex.(i, 3)) => :P 
        )
    else
        tmp =DataFrames.combine(
            groupby(data.loci, [:locus, :population]),
            :genotype => chisq_locus => :chisq
        )
        out_table = select(tmp, :locus, :population, 
            :chisq => (i -> getindex.(i, 1)) => :chisq, 
            :chisq => (i -> getindex.(i, 2)) => :df, 
            :chisq => (i -> getindex.(i, 3)) => :P 
        )
    end
    if correction == "none"
        return out_table
    else
        column_name = Symbol("P_"*correction)
        transform!(out_table, :P => (i -> multitest_missing(i, correction)) => column_name)
    end
end

const hwe = hwe_test
