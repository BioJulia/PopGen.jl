"""
```
ishet(locus::Vector{Union{Missing, NTuple{N,Int8}}} where N)
ishet(locus::Vector{Union{Missing, NTuple{N,Int16}}} where N)
ishet(locus::NTuple{N,Int8} where N)
ishet(locus::NTuple{N,Int16} where N)
ishet(locus::Missing)
```
A series of methods to test if a locus is heterozygous and return `true` if
it is, `false` if it isn't. The vector versions simply broadcast the functions over their
elements.
"""
@inline function ishet(locus::NTuple{N,T} where N where T <: Signed)
    # if allele 1 doesnt equels all others, return true (as in, ishet = true)
    return @inbounds all(locus[1] .== locus)
end

@inline function ishet(locus::Missing)
    # if the locus is Missing, return missing. no muss no fuss
    return missing
end

@inline function ishet(locus::T) where T<:AbstractVector
    ishet.(locus)
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
    Returns the expected heterozygosity of an array of genotypes, calculated
    as 1 - sum of the squared allele frequencies.
"""
function hetero_e(allele_freqs::Vector{T}) where T <: AbstractFloat
    1 - sum(allele_freqs .^ 2)
end



"""
    het_expected(data::PopData)
    Returns the expected heterozygosity of loci in a PopData object.
"""
function het_expected(data::PopData)
    @groupby data.loci :locus {het_exp = hetero_e(allele_freq_vec(:genotype))}
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
        return @groupby data.loci :locus {het_obs = hetero_o(:genotype), het_exp = hetero_e(allele_freq_vec(:genotype))}

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
    het_observed(data::PopData)
Calculate the observed heterozygosity for each locus in `PopData`. Returns a
table of heterozygosity values.
"""
function het_observed(data::PopObj)
    @groupby data.loci :locus {het_obs = hetero_o(:genotype)}
end


"""
    het_population_exp(data::PopObj)
Return a `Dict` of the expected heterozygosity per population for each locus in
a `PopObj`
"""
function het_population_exp(data::PopObj)
    # get population order and store in its own dict key
    popid = unique(data.samples.population |> collect)
    fixed_popid = replace.(popid, " " => "_")
    d = Dict()
    d["pops_exp"] = fixed_popid

    y = deepcopy(data.loci)
    insertcols!(y, 1, :population => data.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        pop_het_vals = Vector{Union{Missing, Float64}}()
        for locus in eachcol(pop[!, :2:end], false)
            a = allele_freq(locus)
            n = length(locus |> skipmissing |> collect)
            het = 1 - (values(a) .^ 2 |> sum)
            het_adjusted = (n/(n-1)) * het
            push!(pop_het_vals, het_adjusted)  # push the sum of hetz to the output array
        end

        # add "_obs" for distinction vs "_exp"
        popname = replace(pop.population[1], " " => "_")*"_exp"
        d[popname] = pop_het_vals
    end
        return d
end

"""
    het_population_obs(data::PopObj)
Return a table of the observed heterozygosity per population for each locus in
`PopData`
"""
function het_population_obs(data::PopObj)
    @groupby data.loci (:locus, :population) {Het_obs = hetero_o(:genotype)}
end


"""
    het_sample(data::PopObj)
Calculate the observed heterozygosity for each individual in `PopData`. Returns
a table of heterozygosity values.
"""
@inline function het_sample(data::PopObj)
    @groupby data.loci (:name) {Het_obs = hetero_o(:genotype)}
end

#TODO
"""
    het_sample(data::PopData, individual::String)
Calculate the observed heterozygosity for an individual in a `PopObj`. Returns
an array of heterozygosity values.
"""
function het_sample(data::PopData, individual::String)
    # calculate observed heterozygosity for an individual
    @where data.loci :name == individual |> @select {het_obs = hetero_o(:genotype)}
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
function hwe_test(data::PopObj; by_pop::Bool = false, correction::String = "none")
    if !by_pop
        output = [
            Vector{Union{Missing,Float64}}(),
            Vector{Union{Missing,Int32}}(),
            Vector{Union{Missing,Float64}}()
            ]
        for locus in eachcol(data.loci, false)
            chisq, df, pval = locus_chi_sq(locus)
            push!.(output, [chisq, df, pval])
        end

        het = heterozygosity(data)
        insertcols!(het, size(het,2)+1, χ² = output[1])
        insertcols!(het, size(het,2)+1, DF = output[2])
        insertcols!(het, size(het,2)+1, P = output[3])

        # corrections
        @info "Χ² test for conformation to Hardy-Weinberg Equilibrium"
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
        popid = unique(data.samples.population |> collect)

        y = deepcopy(data.loci)
        insertcols!(y, 1, :population => data.samples.population)
        y_subdf = groupby(y[!, :], :population) |> collect
        het = heterozygosity(data, "pop")
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
        @info "Χ² test for conformation to Hardy-Weinberg Equilibrium by Population"

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
            p_vals = p_array |> Vector{Union{Missing,Float64}}

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

const hwe = hwe_test


"""
    locus_chi_sq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value
"""
function locus_chi_sq(locus::Vector{<:Union{Missing, NTuple{N,<:Integer}}}) where N
    # Give the function a locus from a genpop object and it will perform the ChiSquared test for HWE
    number_ind = count(i -> i !== missing, locus)

    # Get expected number of genotypes in a locus
    the_allele_dict = allele_freq(locus)

    # split the appropriate pairs into their own vectors
    alleles = Vector{String}()
    p = Vector{Float64}()
    for (i,j) in pairs(the_allele_dict)
        push!(alleles, "$i")
        push!(p, j)
    end

    #Calculate Expected Genotype numbers
    expected_genotype_freq = p * transpose(p) .* number_ind
    expected_genotype_freq = vec(expected_genotype_freq)

    alleles = (alleles .* ",") .* permutedims(alleles)
    alleles = Vector{String}.(sort.(split.(alleles, ",")))
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
    locus_chi_sq(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}},1}) where N
Calculate the chi square statistic and p-value for a locus of a `PopObj` split by
population using `groupby()`. Returns a tuple with chi-square statistic, degrees
of freedom, and p-value.
"""
function locus_chi_sq(locus::SubArray{<:Union{Missing, NTuple{N,<:Integer}},1}) where N
    #Give the function a locus from a genpop object and it will perform the ChiSquared test for HWE
    number_ind = count(i -> i !== missing, locus)

    # split the appropriate pairs into their own vectors
    alleles = Vector{String}()
    p = Vector{Float64}()
    for (i,j) in pairs(the_allele_dict)
        push!(alleles, "$i")
        push!(p, j)
    end

    # create a matrix of allele pair combinations
    alleles = (alleles .* ",") .* permutedims(alleles)
    alleles = Vector{String}.(sort.(split.(alleles, ",")))
    alleles = [parse.(Int16, i) |> Tuple for i in alleles]
    alleles = vec(alleles)

    #Calculate Expected Genotype numbers
    expected_genotype_freq = p * transpose(p) .* number_ind
    expected_genotype_freq = vec(expected_genotype_freq)


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
function multitest_missing(pvals::Array{<:Union{Missing, <:AbstractFloat},1}, correction::String)
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
