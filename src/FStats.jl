function FST(data::PopData; by::String = "global")
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
            HS = mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n))))),
            Het_obs = mean(skipmissing(o)),
            avg_freq = sum(values(avg_allele_freq(alleles)).^2)
            ),
        groupby(het_df, :locus)
    )

    Ht = 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = Ht .- n_df.HS
    DST′ = n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = n_df.HS .+ DST′

    if lowercase(by) != "global"
        return (DST ./ Ht, DST′ ./ HT′)
    else
        Ht = mean(Ht)
        HT′ = mean(HT′)
        DST = mean(DST)
        DST′ = mean(DST′)

        return (DST / Ht, DST′ / HT′,)
    end
end


"""
    f_stat_p(data::PopObj, nperm::Int64)
Performs a permutation test (permuting individuals among populations) to calculate
p-values for F-statistics
"""
function f_stat_p(data::PopData, nperm::Int = 1000)
    n = nperm - 1
    tmp = deepcopy(data)
    observed = FST(data, by = "global")
    perm_dict = Dict{Symbol,Vector{Float64}}(
        :FIS => Vector{Float64}(undef, n),
        :FST => Vector{Float64}(undef, n),
        :FST′ => Vector{Float64}(undef, n)
    )
    for permutation in 1:n
        permute_samples!(tmp)
        tmp_f = FST(tmp, by = "global")
        perm_dict[:FIS][permutation] = tmp_f[:FIS]
        perm_dict[:FST][permutation] = tmp_f[:FST]
        perm_dict[:FST′][permutation] = tmp_f[:FST′]
    end

    for (k,v) in perm_dict
        observed[k] = (1 + sum(observed[k] .<= v))/nperm
    end

    return observed
    #return map(i -> sum((observed .<= sims)[!,i])/nperm, 1:ncol(observed))
end
