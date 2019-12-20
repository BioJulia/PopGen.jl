"""
    Fstats(x::PopObj, mode::String = "global")
Returns a DataFrame of the F-statistics of a `PopObj`. See
`FST` to perform pairwise FST comparisons between populations.

### mode options
- `global` : provides FIT, FIS, and FST for an entire population (default)
- `loci` : provides FIS, FIT, and FST for each locus
- `sample` : provides FIS and FIT for each sample
"""
function Fstats(x::PopObj, mode::String = "global")
    hets = He(x, "pop")
    H_t = het_expected(x) |> mean
    H_s = map(i -> mean(i.het_exp), hets)
    H_i = map(i -> mean(i.het_obs), hets)

    if mode == "locus" || mode == "loci"
        @info "F-statistics by Locus"

        # FIT = (HT- HI)/HT
        FIT = map(i -> (H_t - i) / H_t, het_observed(x))

        # FIS = (HS - HI) / HS
        HI = het_observed(x)
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
        for i = 1:length(loci(x))
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

    elseif mode == "sample" || mode == "ind" || mode == "individual"
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
        FIT = map(i -> (H_t - i) / H_t, H_i_sample.het)

        return DataFrame(
            name = samples(x),
            FIS = FIS,
            FIT = round.(FIT, digits = 4),
        )

    elseif mode == "global"
        H_i_sample = het_sample(x)
        @info "Global F-statistics"

        # FIS = (HS - HI) / HS
        FIS = (H_s .- H_i) ./ H_s |> mean

        # FIT = (HT- HI)/HT
        FIT = map(i -> (H_t - i) / H_t, H_i_sample.het) |> mean

        # FST = (HT- HS) / HT
        FST = map(i -> (H_t - i) / H_t, H_s) |> mean

        return DataFrame(
            FIS = round.(FIS, digits = 4),
            FIT = round.(FIT, digits = 4),
            FST = round.(FST, digits = 4),
        )
    else
        error("please specify \"global\", \"loci\", or \"sample\" as the second argument")
    end
end

function fst_global(data:: PopObj)
    H_i_sample = het_sample(data)
    @info "Global F-statistics"

    # FIS = (HS - HI) / HS
    FIS = (H_s .- H_i) ./ H_s |> mean

    # FIT = (HT- HI)/HT
    FIT = map(i -> (H_t - i) / H_t, H_i_sample.het) |> mean

    # FST = (HT- HS) / HT
    FST = map(i -> (H_t - i) / H_t, H_s) |> mean

    return DataFrame(
        FIS = round.(FIS, digits = 4),
        FIT = round.(FIT, digits = 4),
        FST = round.(FST, digits = 4),
    )
end



#======= FST_alpha =======#
hets = He(x, "pop")
H_s = map(i -> mean(i.het_exp), hets)
#H_i = map(i -> mean(i.het_obs), hets)
