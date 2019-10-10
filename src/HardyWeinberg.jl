"""
    het_observed(x::PopObj)
Calculate the observed heterozygosity for each locus in a `PopObj`. Returns an
array of heterozygosity values.
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
Calculate the expected heterozygosity for each locus in a `PopObj`. Returns an
array of heterozygosity values.
"""
function het_expected(x::PopObj)
    het_vals = []
    for locus in eachcol(x.loci, false)
        a = allele_freq(locus) # get allele freqs at locus
        a = a |> values |> collect # isolate the freqs
        homz = a .^ 2 |> sum
        hetz = 1 - homz
        push!(het_vals, hetz)  # push the sum of hetz to the output array
    end
    return het_vals |> Array{Float64,1}
end

"""
    het_sample(x::PopObj)
Calculate the observed heterozygosity for each individual in a `PopObj`. Returns
an array of heterozygosity values.
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
Return a Dict of the observed heterozygosity per population for each locus in
a `PopObj`
"""
function het_population_obs(x::PopObj)
    # get population order and store in its own dict key
    popid = unique(x.samples.population |> collect)
    fixed_popid = [replace(i, " " => "") for i in popid]
    d = Dict()
    d["pops_obs"] = fixed_popid

    y = deepcopy(x.loci)
    insertcols!(y, 1, :population => x.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        pop_het_vals = []
        for locus in eachcol(pop[!, :2:end], false)
            # get genotype freqs at locus
            a = geno_freq(locus)

            # add condition if the locus is missing from that subpop
            if a === missing
                push!(pop_het_vals, missing)
                continue
            end
            tmp = 0
            for geno in collect(keys(a))
                # create hom geno
                geno_hom = fill(geno[1], length(geno)) |> Tuple
                if geno != geno_hom    # test if geno isn't homozygous
                    tmp += a[geno]     # if true, add freq to total in tmp
                end
            end
            push!(pop_het_vals, tmp)
        end

        # add "_obs" for distinction vs "_exp"
        popname = replace(pop.population[1], " " => "")*"_obs"
        d[popname] = pop_het_vals |> Array{Union{Missing,Float64},1}
    end
    return d
end


"""
    het_population_exp(x::PopObj)
Return a `Dict` of the expected heterozygosity per population for each locus in
a `PopObj`
"""
function het_population_exp(x::PopObj)
    # get population order and store in its own dict key
    popid = unique(x.samples.population |> collect)
    fixed_popid = [replace(i, " " => "") for i in popid]
    d = Dict()
    d["pops_exp"] = fixed_popid

    y = deepcopy(x.loci)
    insertcols!(y, 1, :population => x.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        pop_het_vals = []
        for locus in eachcol(pop[!, :2:end], false)
            a = allele_freq(locus) # get allele freqs at locus
            a = a |> values |> collect # isolate the freqs
            homz = a .^ 2 |> sum
            hetz = 1 - homz
            push!(pop_het_vals, hetz)  # push the sum of hetz to the output array
        end

        # add "_obs" for distinction vs "_exp"
        popname = replace(pop.population[1], " " => "")*"_exp"
        d[popname] = pop_het_vals |> Array{Union{Missing,Float64},1}
    end
        return d
end


"""
    heterozygosity(x::PopObj, mode = "locus")
Calculate observed and expected heterozygosity in a `PopObj`
## Example
heterozygosity(nancycats(), "population" )

### Modes
- `"locus"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
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
        obs = het_population_obs(x)
        exp = het_population_exp(x)
        merged_het = merge(obs, exp)
        locinames = String.(names(x.loci))

        # remake obs/exp distinction
        obs_pops = obs["pops_obs"] .* "_obs"
        exp_pops = exp["pops_exp"] .* "_exp"

        # interleave pop_obs/exp for column names
        #col_names = Iterators.flatten(zip(obs_pops, exp_pops)) |> collect
        df = DataFrame(locus = locinames)
        out_df = []
        for (i,j) in zip(obs_pops, exp_pops)
            tmp = DataFrame(locus = locinames, het_obs = merged_het[i], het_exp = merged_het[j])
            push!(out_df, tmp)
        end
        #for pop_het in col_names
        #    insertcols!(df, size(df,2)+1, Symbol(pop_het) => merged_het[pop_het])
        #end
    end
    return Tuple(out_df)
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
    the_allele_dict = allele_freq(locus)
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
    locus_chi_sq(locus::SubArray{Union{Missing, Tuple},1})
Calculate the chi square statistic and p-value for a locus of a `PopObj` split by
population using `groupby()`. Returns a tuple with chi-square statistic, degrees
of freedom, and p-value.
"""
function locus_chi_sq(locus::SubArray{Union{Missing,Tuple},1})
    #Give the function a locus from a genpop object and it will perform the ChiSquared test for HWE
    number_ind = count(i -> i !== missing, locus)

    ## Get expected number of genotypes in a locus
    the_allele_dict = allele_freq(locus)
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
    hwe_test(x::PopObj; by_pop::Bool = false; correction = "none")
Calculate chi-squared test of HWE for each locus and returns observed and
expected heterozygosity with chi-squared, degrees of freedom and p-values
for each locus. Use `by_pop = true` to perform this separately for each
population (default: by_pop = false). Use `correction =` to specify a P-value
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
function hwe_test(x::PopObj; by_pop::Bool = false, correction::String = "none")
    if by_pop == false
        output = [[], [], []]
        for locus in eachcol(x.loci, false)
            chisq, df, pval = locus_chi_sq(locus)
            push!.(output, [chisq, df, pval])
        end
        #output = map(locus_chi_sq, eachcol(x.loci, false)) |> DataFrame
        het = heterozygosity(x)
        insertcols!(het, size(het,2)+1, χ² = output[1] |> Array{Union{Missing,Float64},1})
        insertcols!(het, size(het,2)+1, DF = output[2] |> Array{Union{Missing,Float64},1})
        insertcols!(het, size(het,2)+1, P = output[3] |> Array{Union{Missing,Float64},1})
        ## corrections
        if correction == "none"
            return het
        else
            corrected_P = multitest_missing(het.P, correction)

            # add the adjusted p-vals back to the dataframe and return
            insertcols!(het, 7, Pcorr = corrected_P )
        end
    else
        out_array = []
        # get population order and store in its own dict key
        popid = unique(x.samples.population |> collect)

        y = deepcopy(x.loci)
        insertcols!(y, 1, :population => x.samples.population)
        y_subdf = groupby(y[!, :], :population) |> collect
        het = heterozygosity(x, "pop")
        for (eachpop, het_pop) in zip(y_subdf, het)
            output = [[], [], []]
            for locus in eachcol(eachpop, false)[2:end]
                chisq, df, pval = locus_chi_sq(locus)
                push!.(output, [chisq, df, pval])
            end
            insertcols!(het_pop, size(het_pop,2)+1, χ² = output[1] |> Array{Union{Missing,Float64},1})
            insertcols!(het_pop, size(het_pop,2)+1, DF = output[2] |> Array{Union{Missing,Float64},1})
            insertcols!(het_pop, size(het_pop,2)+1, P = output[3] |> Array{Union{Missing,Float64},1})
            push!(out_array, het_pop)
        end

        ## corrections
        if correction == "none"
            d = Dict([Symbol(i) => j  for (i,j) in zip(popid, out_array)])
            return NamedTuple{Tuple(Symbol.(popid))}(values(d))
        else
            p_array = []

            #for each in out_array
            for i in out_array
                append!(p_array, i.P)
            end

            # convert p_array to suitible Type
            p_vals = p_array |> Array{Union{Missing,Float64},1}

            # do actual correction
            corrected_P = multitest_missing(p_vals, correction)

            # split corrected pvals by #loci
            p_corrs = Iterators.partition(corrected_P, size(out_array[1],1)) |> collect

            # add each segment of split-pvals to the appropriate df
            for i in 1:length(out_array)
                insertcols!(out_array[i], size(out_array[i],2)+1, Pcorr = p_corrs[i])
            end
            d = Dict([Symbol(i) => j  for (i,j) in zip(popid, out_array)])
            return NamedTuple{Tuple(Symbol.(popid))}(values(d))
        end
    end
end


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
function multitest_missing(pvals::Array{Union{Missing, Float64},1}, correction::String)
    # make seperate array for non-missing P vals
    p_no_miss = skipmissing(pvals) |> collect
    # get indices of where original missing are
    miss_idx = findall(i -> i === missing, pvals)

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

    correct = adjust(p_no_miss, d[lowercase(correction)]) |> Array{Union{Missing, Float64},1}

    # re-add missing to original positions
    for i in miss_idx
        insert!(correct, i, missing)
    end
    return correct
end
