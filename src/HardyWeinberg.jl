"""
    het_observed(x::PopObj)
Calculate the observed heterozygosity for each locus in a `PopObj`
"""
function het_observed(x::PopObj)
    het_vals = []
    for locus in eachcol(x.loci, false)
        a = geno_freq(locus)  # get genotype freqs at locus
        tmp = 0
        for geno in collect(keys(a))
            geno_hom = fill(geno[1], length(geno)) |> Tuple   # create hom geno
            if geno != geno_hom        # test if geno isn't homozygous
                tmp += a[geno]     # if true, add freq to total in tmp
            end
        end
        push!(het_vals, tmp)
    end
    return (het_vals) |> Array{Float64,1}
end

"""
    het_expected(x::PopObj)
Calculate the expected heterozygosity for each locus in a `PopObj`
"""
function het_expected(x::PopObj)
    het_vals = []
    for locus in eachcol(x.loci, false)
        a = allele_freq_mini(locus) # get allele freqs at locus
        a = a |> values |> collect # isolate the freqs
        homz = a .^ 2 |> sum
        hetz = 1 - homz
        push!(het_vals, hetz)  # push the sum of hetz to the output array
    end
    return het_vals |> Array{Float64,1}
end

"""
    het_sample(x::PopObj)
Calculate the observed heterozygosity for each individual in a `PopObj`
"""
function het_sample(x::PopObj)
    #transpose loci df to have samples as columns
    y = deepcopy(x.loci)
    insertcols!(y, 1, :name => x.samples.name)
    insertcols!(y, 2, :id => 1:length(y[!, 1]))
    tmp_stack = stack(y, names(y[!, 3:end]))
    by_ind = unstack(tmp_stack, :variable, :id, :value)
    select!(by_ind, Not(:variable))
    # calculate observed heterozygosity like for loci
    het_vals = []
    for samp in eachcol(by_ind, false)
        a = geno_freq(samp)  # get genotype freqs at sample
        #delete!(a, missing)     # remove missing values
        tmp = 0
        for geno in collect(keys(a))
            geno_hom = fill(geno[1], length(geno)) |> Tuple   # create hom geno
            if geno != geno_hom        # test if geno isn't homozygous
                tmp += a[geno]     # if true, add freq to total in tmp
            end
        end
        push!(het_vals, tmp)
    end
    return DataFrame(
        name = x.samples.name,
        het = (het_vals |> Array{Float64,1}),
    )
end

"""
    het_population_obs(x::PopObj)
Return the observed heterozygosity per population for each locus in a `PopObj`
"""
function het_population_obs(x::PopObj)
    d = Dict()
    popnames = []
    y = deepcopy(x.loci)
    insertcols!(y, 1, :population => x.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        pop_het_vals = []
        for locus in eachcol(pop[!, :2:end], false)
            a = geno_freq(locus)  # get genotype freqs at locus
            # add condition if the locus is missing from that subpop
            if a === missing
                push!(pop_het_vals, missing)
                continue
            end
            tmp = 0
            for geno in collect(keys(a))
                geno_hom = fill(geno[1], length(geno)) |> Tuple   # create hom geno
                if geno != geno_hom        # test if geno isn't homozygous
                    tmp += a[geno]     # if true, add freq to total in tmp
                end
            end
            push!(pop_het_vals, tmp)
        end
        # convert to include missing
        pop_het_conv = pop_het_vals |> Array{Union{Missing,Float64},1}
        # get the population name and remove whitespaces
        popname = replace(pop.population[1], " " => "")
        d[popname] = pop_het_conv
        push!(popnames, popname)
    end
    #return_df = DataFrame(d)[!,Symbol.(popnames)]
    insertcols!(DataFrame(d), 1, :locus => string.(x.loci |> names))
end


"""
    het_population_exp(x::PopObj)
Return the expected heterozygosity per population for each locus in a `PopObj`
"""
function het_population_exp(x::PopObj)
    d = Dict()
    popnames = []
    y = deepcopy(x.loci)
    insertcols!(y, 1, :population => x.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        pop_het_vals = []
        for locus in eachcol(pop[!, :2:end], false)
            a = allele_freq_mini(locus) # get allele freqs at locus
            a = a |> values |> collect # isolate the freqs
            homz = a .^ 2 |> sum
            hetz = 1 - homz
            push!(pop_het_vals, hetz)  # push the sum of hetz to the output array
        end
        #het_vals |> Array{Float64,1}
        # get the population name and remove whitespaces
        popname = replace(pop.population[1], " " => "")
        d[popname] = pop_het_vals
        push!(popnames, popname)
    end
    #return_df = DataFrame(d)[!,Symbol.(popnames)]
    insertcols!(DataFrame(d), 1, :locus => string.(x.loci |> names))
end


"""
    heterozygosity(x::PopObj, mode = "locus")
Calculate observed and expected heterozygosity in a `PopObj`
## Example
heterozygosity(nancycats, "population" )

### Modes
- `"locus"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual
- `"population"` or `"pop"` : heterozygosity per population (PopObj.samples.population)
"""
function heterozygosity(x::PopObj, mode = "locus")
    if mode == "locus"
        obs = het_observed(x)
        exp = het_expected(x)
        locinames = String.(names(x.loci))
        return DataFrame(locus = locinames, het_obs = obs, het_exp = exp)
    elseif lowercase(mode) == "sample" || lowercase(mode) == "ind" || lowercase(mode) == "individual"
        return het_sample(x)
    elseif lowercase(mode) == "pop" || lowercase(mode) == "population"
        return het_population_obs(x)
    end
end

const het = heterozygosity
const He = heterozygosity

"""
    locus_chi_sq(locus::Array{Union{Missing, Tuple},1})
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value
"""
function locus_chi_sq(locus::Array{Union{Missing,Tuple},1})
    #Give the function a locus from a genpop object and it will perform the ChiSquared test for HWE
    number_ind = count(i -> i !== missing, locus)

    ## Get expected number of genotypes in a locus
    the_allele_dict = allele_freq_mini(locus)
    p = the_allele_dict |> values |> collect

    #Calculate Expected Genotype numbers
    expected_genotype_freq = p * transpose(p) .* number_ind
    expected_genotype_freq = vec(expected_genotype_freq)

    alleles = ["$i" for i in the_allele_dict |> keys |> collect]
    alleles = (alleles .* ",") .* permutedims(alleles)
    alleles = Array{String}.(sort.(split.(alleles, ",")))
    alleles = [parse.(Int16, i) |> Tuple for i in alleles]
    alleles = vec(alleles)

    expected = Dict()
    for (geno, freq) in zip(alleles, expected_genotype_freq)
        expected[geno] = get!(expected, geno, 0) + freq
    end

    ## Get observed number of genotypes in a locus
    observed = geno_freq(locus)
    for j in keys(observed)
        observed[j] = observed[j] * number_ind
    end

    chisq_stat = expected
    for genotype in keys(expected)
        o = get(observed, genotype, 0)
        e = get(expected, genotype, 0)

        chisq_stat[genotype] = (o - e)^2 / e
    end
    chisq_stat = values(chisq_stat) |> sum
    df = (length(alleles) - length(the_allele_dict)) / 2

    if df > 0
        chi_sq_dist = Distributions.Chisq(df)
        p_val = 1 - cdf(chi_sq_dist, chisq_stat)
    else
        p_val = missing
    end
    return (chisq_stat, df, p_val)
end


"""
    hwe_test(x::PopObj; correction = "none")
Calculate chi-squared test of HWE for each locus and returns observed and
exected heterozygosity with chi-squared, degrees of freedom and p-values
for each locus.  Use `correction =` to specify a P-value correction method
for multiple testing.

#### example
`hwe_test(gulfsharks(), correction = "bh")`

### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"benjamini"` : Benjamini adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` or `"b-h"` : Benjamini-Hochberg adjustment
- `"by"` or `"b-y"`: Benjamini-Yekutieli adjustment
- `"bl"` or `"b-l"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forward stop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` or `"b-c"` : Barber-Candès adjustment
"""
function hwe_test(x::PopObj; correction = "none")
    output = [[], [], []]
    for locus in eachcol(x.loci, false)
        chisq, df, pval = locus_chi_sq(locus)
        push!.(output, [chisq, df, pval])
    end
    #output = map(locus_chi_sq, eachcol(x.loci, false)) |> DataFrame
    het = heterozygosity(x)
    insertcols!(het, 4, ChiSq = output[1] |> Array{Union{Missing,Float64},1})
    insertcols!(het, 5, DF = output[2] |> Array{Union{Missing,Float64},1})
    insertcols!(het, 6, P = output[3] |> Array{Union{Missing,Float64},1})
    ## corrections
    if correction == "none"
        return het
    else
        # make seperate array for non-missing P vals
        p_no_miss = skipmissing(het.P) |> collect
        # get indices of where original missing are
        miss_idx = findall(i -> i === missing, het.P)

        # make a dict of all possible tests and their respective functions
        d = Dict(
            "bonferroni" => Bonferroni(),
            "holm" => Holm(),
            "benjamini" => Benjamini(),
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

        correct = adjust(p_no_miss, d[lowercase(correction)]) |> Array{Any,1}

        # re-add missing to original positions
        for i in miss_idx
            insert!(correct, i, missing)
        end

        # add the adjusted p-vals back to the dataframe and return
        insertcols!(het, 7, Pcorr = correct |> Array{Union{Missing,Float64},1})
    end
end
