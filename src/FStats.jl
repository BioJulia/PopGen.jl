"""
    f_stats(x::PopObj, mode::String = "global")
Returns a DataFrame of the F-statistics of a `PopObj`. See
`FST` to perform pairwise FST comparisons between populations.

### mode options
- `global` : provides FIT, FIS, and FST for an entire population (default)
- `loci` : provides FIS, FIT, and FST for each locus
- `sample` : provides FIS and FIT for each sample
"""
function f_stats(x::PopObj, mode::String = "global")

    if mode == "locus" || mode == "loci"

        fst_locus(x)

    elseif mode == "sample" || mode == "ind" || mode == "individual"

        fst_sample(x)

    elseif mode == "global"

        fst_global(x)

    else
        error("please specify \"global\", \"loci\", or \"sample\" as the second argument")
    end
end


"""
    fst_global(data::PopObj, round::Bool = true)
Calculates the three common F-statistics Fᵢₛ, Fᵢₜ, and Fₛₜ on a PopObj, adjusted for
different sample sizes in populations. If `round` is `true` (default), the results are
rounded to 4 decimal places, otherwise `false` returns the entire values as calculated.
"""
function fst_global(data::PopObj, round::Bool = true)

    ## H_i = (Hobs1*N1 + Hobs2*N2 + etc..) / Ntotal
    popcounts = populations(data).count
    pop_hets = het(data, "pop")
    obs = map(i -> mean(i.het_obs), pop_hets) |> collect
    exp = map(i -> mean(i.het_exp), pop_hets) |> collect
    H_i = sum(obs .* popcounts) / sum(popcounts)

    ## n/n-1 correction
    correction = sum(popcounts) / (sum(popcounts) - 1)
    ## H_s = (Hexp1*N1 + Hexp2*N2 + etc..) / Ntotal
    H_s = (sum(exp .* popcounts) / sum(popcounts)) * correction

    ## H_t = 1 - ∑(p²+q²)     (expected Het)
    H_t = mean(het(data).het_exp)

    @info "Global F-statistics"

    # FIS = (HS - HI) / HS
    FIS = (H_s - H_i) / H_s

    # FIT = (HT- HI)/HT
    FIT = (H_t - H_i) / H_t

    # FST = (HT- HS) / HT
    FST = (H_t - H_s) / H_t

    if round == false
        return DataFrame(
            FIS = FIS,
            FIT = FIT,
            FST = FST
        )
    else
        return DataFrame(
            FIS = Base.round.(FIS, digits = 4),
            FIT = Base.round.(FIT, digits = 4),
            FST = Base.round.(FST, digits = 4),
        )
    end
end





function f_pop(data::PopData)
    hets = heterozygosity(data, "population")
    transform(hets, F => (:het_exp - :het_obs)/:het_exp)
end




function fst_locus(data::PopObj)
    hets = He(data, "pop")
    H_t = het_expected(data) |> mean
    H_s = map(i -> mean(i.het_exp), hets)
    H_i = map(i -> mean(i.het_obs), hets)

    @info "F-statistics by Locus"

    # FIT = (HT- HI)/HT
    FIT = map(i -> (H_t - i) / H_t, het_observed(data))

    # FIS = (HS - HI) / HS
    HI = het_observed(data)
    HS = map(i -> i.het_exp, hets)
    FIS = []
    for (i, j) in enumerate(HI)
        locus = map(a -> a[i], HS)
        fis = (locus .- j) ./ locus |> mean
        push!(FIS, fis)
    end
    FIS = round.(FIS, digits = 4) |> Array{Float64,1}

    # FST = (HT- HS) / HT
    FST = []
    for i = 1:length(loci(data))
        locus = map(c -> c[i], HS)
        fst = (H_t .- locus) ./ H_t |> mean
        push!(FST, fst)
    end
    FST = round.(FST, digits = 4) |> Array{Float64,1}

    return DataFrame(
        locus = loci(x),
        FIS = FIS,
        FIT = round.(FIT, digits = 4),
        FST = FST,
    )
end


function fst_sample(data::PopObj)
    H_i_sample = het_sample(x)
    @info "F-statistics by Sample"
    # split sample heterozygosities into their own GroupedDataFrames
    insertcols!(H_i_sample, 2, :population => x.samples.population)
    H_i_sample_pop = groupby(H_i_sample[!, :], :population)

    FIS = []
    for (i, j) in enumerate(H_i_sample_pop)
        # FIS = (HS - HI) / HS
        tmp = map(k -> (H_s[i] - k) / H_s[i], j.het)
        append!(FIS, tmp)
    end

    FIS = round.(FIS, digits = 4) |> Array{Float64,1}

    # FIT = (HT- HI)/HT
    # H_t =
    FIT = map(i -> (H_t - i) / H_t, H_i_sample.het)

    return DataFrame(
        name = samples(x),
        FIS = FIS,
        FIT = round.(FIT, digits = 4),
    )
end

"""
    f_stat_p(data::PopObj, nperm::Int64)
Performs a permutation test (permuting individuals among populations) to calculate
p-values for F-statistics
"""
function f_stat_p(data::PopObj, nperm::Int64 = 1000)
    tmp = deepcopy(data)

    observed = f_stats(data, "global")
    sims = deepcopy(observed)
    for boot in 1:(nperm-1)
        shuffle!(tmp.samples.population)
        append!(sims, f_stats(tmp, "global"))
    end
    return map(i -> sum((observed .<= sims)[!,i])/nperm, 1:ncol(observed))
end
