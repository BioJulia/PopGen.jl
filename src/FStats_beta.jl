function fstat(data::PopData; by::String = "global")
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data.loci, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )
    # collapse down to retrieve averages and counts
    n_df = DataFrames.combine(
        [:het_obs, :het_exp, :n, :alleles] => (o,e,n,alleles) -> (
            count = sum(.!iszero.(n)),
            mn = sum(.!iszero.(n)) ./ sum(reciprocal.(n)),
            Hs_nei = mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n))))),
            Het_obs = mean(skipmissing(o)),
            avg_freq = sum(values(avg_allele_freq(alleles)).^2)
            ),
        groupby(het_df, :locus)
    )

    Ht = 1.0 .- n_df.avg_freq .+ (n_df.Hs_nei ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    FIS =  1.0 .- (n_df.Het_obs ./ n_df.Hs_nei)
    DST = Ht .- n_df.Hs_nei
    FST = DST ./ Ht
    DSTP = n_df.count ./ (n_df.count .- 1) .* DST
    DEST = DSTP ./ (1 .- n_df.Hs_nei)
    HTP = n_df.Hs_nei .+ DSTP
    FSTP = DSTP ./ Ht

    if lowercase(by) == "global"
        Het_obs = mean(n_df.Het_obs)
        Hs_nei = mean(n_df.Hs_nei)
        Ht = mean(Ht)
        FIS = mean(FIS)
        DST = mean(DST)
        FST = mean(FST)
        DSTP = mean(DSTP)
        DEST = mean(DEST)
        HTP = mean(HTP)
        FSTP = mean(FSTP)

        DataFrame(
        :Het_obs => round.(Het_obs, digits = 4),
        :Hs_nei => round.(Hs_nei, digits = 4),
        :HT => round.(Ht , digits = 4),
        :DST => round.(DST, digits = 4),
        :HTP => round.(HTP, digits = 4),
        :DSTP => round.(DSTP, digits = 4),
        :FST => round(DST/Ht, digits = 4),
        :FSTP => round(DSTP/HTP, digits = 4),
        :FIS => round(1 - (Het_obs / Hs_nei), digits = 4),
        :DEST => round(DSTP / (1 - Hs_nei), digits = 4)
        )
    elseif by ∈ ["locus", "loci"]
        insertcols!(
        select!(
        n_df, :locus,
        :Het_obs => (i -> round.(i, digits = 4)) => :Het_obs,
        :Hs_nei => (i -> round.(i, digits = 4)) => :Hs_nei,
        ),
        :HT => round.(Ht, digits = 4),
        :DST => round.(DST, digits = 4),
        :HTP => round.(HTP, digits = 4),
        :DSTP => round.(DSTP, digits = 4),
        :FST => round.(FST, digits = 4),
        :FSTP => round.(FSTP, digits = 4),
        :FIS => round.(FIS, digits = 4),
        :DEST => round.(DEST, digits = 4)
        )
    end
end
