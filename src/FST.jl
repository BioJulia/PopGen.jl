function Fstats(x::PopObj, mode::String = "global")
    hets = He(x, "pop")
    H_t = het_expected(x) |> mean
    H_s = map(i -> mean(i.het_exp), hets)
    H_i = map(i -> mean(i.het_obs), hets)
    if mode == "locus" || mode == "loci"
        @info "F-statistics by Locus"

        # FIT = the locus observed mean compared to the total expected
        # FIT = (HT- HI)/HT
        FIT = map(i -> (H_t - i) / H_t, het_observed(x))

        # the locus expected vs total
        # FST = (HT- HS) / HT
        FST = map(i -> (H_t - i) / H_t, H_s)
        println("length of FST = ", length(FST))

        return DataFrame(locus = loci(x),
                #FIS = FIS,
                FIT = round.(FIT, digits = 4),
                FST = round.(FST, digits = 4)
                )

    elseif mode == "sample" || mode == "ind" || mode == "individual"
        @info "F-statistics by Sample"
        H_i_sample = het_sample(x)
        # split sample heterozygosities into their own GroupedDataFrames
        insertcols!(H_i_sample, 2, :population => x.samples.population)
        H_i_sample_pop = groupby(H_i_sample[!, :], :population)

        FIS = []
        for (i,j) in enumerate(H_i_sample_pop)
            # FIS = (HS - HI) / HS
            tmp = map(k -> (H_s[i] - k) / H_s[i], j.het)
            append!(FIS, tmp)
        end

        FIS = round.(FIS, digits = 4) |> Array{Float64,1}

        # FIT = (HT- HI)/HT
        FIT = map(i -> (H_t - i) / H_t, H_i_sample.het)
        return DataFrame(name = samples(x),
                FIS = FIS,
                FIT = round.(FIT, digits = 4)
                )

    elseif mode == "global"
        @info "Global F-statistics"

        # FIS = (HS - HI) / HS
        FIS = (H_s .- H_i) ./ H_s |> mean

        # FIT = (HT- HI)/HT
        FIT = map(i -> (H_t - i) / H_t, H_i_sample.het) |> mean

        # FST = (HT- HS) / HT
        FST = map(i -> (H_t - i) / H_t, H_s) |> mean

        return DataFrame(FIS = round.(FIS, digits = 4),
                    FIT = round.(FIT, digits = 4),
                    FST = round.(FST, digits = 4)
                )
    end
end
