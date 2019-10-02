"""
    geno_freq_alpha(x::Array{Union{Missing, Tuple},1})
Calculate genotype frequencies of all loci in a `PopObj`
"""
function geno_freq_alpha(x::Array{Union{Missing,Tuple},1})
    d = Dict()
    for row in x
        # sum up missing
        if row === missing
            continue
        else
        # sum up non-missing genotypes
            d[row] = get!(d, row, 0) + 1
        end
    end
    total = values(d) |> sum    # sum of all non-missing genotypes
    [d[i] = d[i] / total for i in keys(d)] # genotype count/total
    return d
end


"""
  allele_freq_beta(x::PopObj)
Returns an array of `Dicts` of allele counts per locus
"""
function allele_freq_alpha(x::PopObj, include_missing = false)
    y = PopOpt(x)
    tmp = names(y.loci)[1]  # restrict to single locus for testing
    d = Dict()
    for row in y.loci[!, tmp]
        if row === missing
            continue
        else
            for allele in row
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum
    [d[i] = d[i] / total for i in keys(d)]
    return d
end


"""
    allele_freq_mini(x::Array{Union{Missing, Tuple},1})
Calculate allele counts for a single locus of a `PopObj`
"""
function allele_freq_mini(x::Array{Union{Missing,Tuple},1})
    d = Dict()
    for row in x
        # sum up missing
        if row === missing
            continue
        else
        # sum up alleles
            for allele in row
                d[allele] = get!(d, allele, 0) + 1
            end
        end
    end
    total = values(d) |> sum    # sum of all non-missing alleles
    [d[i] = d[i] / total for i in keys(d)]  # allele count / sum
    return d
end


"""
    het_observed(x::PopObj)
Calculate the observed heterozygosity for each locus in a `PopObj`
"""
function het_observed(x::PopObj)
    het_vals = []
    for locus in eachcol(x.loci, false)
        a = geno_freq_alpha(locus)  # get genotype freqs at locus
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
    return (het_vals) |> Array{Float64,1}
    #locinames = String.(names(x.loci))
    #return DataFrame(locus = locinames, het_obs = Array{Float64,1}(het_vals))
end

"""
    het_expected(x::PopObj)
Calculate the expected heterozygosity for each locus in a `PopObj`
"""
function het_expected(x::PopObj)
    het_vals = []
    for locus in eachcol(x.loci, false)
        a = allele_freq_mini(locus) # get allele freqs at locus
        #delete!(a, missing)     # remove missing values
        a = a |> values |> collect # isolate the freqs
        homz = a .^ 2 |> sum
        hetz = 1 - homz
        push!(het_vals, hetz)  # push the sum of hetz to the output array
    end
    return het_vals |> Array{Float64,1}
    #locinames = String.(names(x.loci))
    #return DataFrame(locus = locinames, het_exp = Array{Float64,1}(het_vals))
end

"""
    het_sample(x::PopObj)
Calculate the heterozygosity for each individual in a `PopObj`
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
        a = geno_freq_alpha(samp)  # get genotype freqs at sample
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

function het_population_obs(x::PopObj)
    pop_het_vals = []
    y = deepcopy(x)
    insertcols!(y, 1, :population => x.samples.population)
    y_subdf = groupby(y[!, :], :population)
    for pop in y_subdf
        for locus in eachcol(pop, false)
            a = geno_freq_alpha(locus)  # get genotype freqs at locus
            tmp = 0
            for geno in collect(keys(a))
                geno_hom = fill(geno[1], length(geno)) |> Tuple   # create hom geno
                if geno != geno_hom        # test if geno isn't homozygous
                    tmp += a[geno]     # if true, add freq to total in tmp
                end
            end
            push!(het_vals, tmp)
        end
    end
    return (het_vals) |> Array{Float64,1}
end

"""
    heterozygosity(x::PopObj)
Calculate observed and expected heterozygosity in a `PopObj`.

### Modes
- `"locus"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` : heterozygosity per individual
- `"population"` or `"pop"` : heterozygosity per population (PopObj.samples.population)
"""
function heterozygosity(x::PopObj, mode = "locus")
    if mode == "locus"
        obs = het_observed(x)
        exp = het_expected(x)
        locinames = String.(names(x.loci))
        return DataFrame(locus = locinames, het_obs = obs, het_exp = exp)
    elseif mode == "sample" || "ind" || "individual"
        return het_sample(x)
    elseif mode == "pop" || "population"
        return
    end
end

const het = heterozygosity
const He = heterozygosity
