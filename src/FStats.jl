function FST_global(data::PopData)
    # observed/expected het per locus per pop
    het_df = @inbounds DataFrames.combine(
        groupby(data.loci, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )

    # collapse down to retrieve averages and counts
    n_df = @inbounds DataFrames.combine(
        groupby(het_df, :locus),
        :n => (n -> sum(.!iszero.(n))) => :count,
        :n => (n -> sum(.!iszero.(n)) ./ sum(reciprocal.(n))) => :mn,
        [:het_obs, :het_exp, :n] => ((o, e, n) -> mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n)))))) => :HS,
        :het_obs => (o -> mean(skipmissing(o))) => :Het_obs,
        :alleles => (alleles -> sum(values(avg_allele_freq(alleles)).^2)) => :avg_freq
    )

    HT = mean(@inbounds 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count))
    DST = mean(@inbounds HT .- n_df.HS)
    DST′ = mean(@inbounds n_df.count ./ (n_df.count .- 1) .* DST)
    HT′ = mean(@inbounds n_df.HS .+ DST′)

    return Dict{Symbol, Float64}(
        :FST => DST / HT,
        :FST′ => DST′ / HT′
    )
end

function FST_locus(data::PopData)
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

    HT = 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = HT .- n_df.HS
    DST′ = n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = n_df.HS .+ DST′

    return Dict{Symbol, Vector{Float64}}(
        :FST => DST ./ HT,
        :FST′ => DST′ ./ HT′
        )
end

function FST_global(data::AbstractDataFrame)
    # observed/expected het per locus per pop
    het_df = @inbounds DataFrames.combine(
        groupby(data, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => (i -> hetero_e(i)) => :het_exp,
        :genotype => allele_freq => :alleles
    )

    # collapse down to retrieve averages and counts
    n_df = @inbounds DataFrames.combine(
        groupby(het_df, :locus),
        :n => (n -> sum(.!iszero.(n))) => :count,
        :n => (n -> sum(.!iszero.(n)) ./ sum(reciprocal.(n))) => :mn,
        [:het_obs, :het_exp, :n] => ((o, e, n) -> mean(skipmissing(gene_diversity_nei87.(e,o,sum(.!iszero.(n)) ./ sum(reciprocal.(n)))))) => :HS,
        :het_obs => (o -> mean(skipmissing(o))) => :Het_obs,
        :alleles => (alleles -> sum(values(avg_allele_freq(alleles)).^2)) => :avg_freq
    )


    HT = mean(@inbounds 1.0 .- n_df.avg_freq + (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count))
    DST = mean(@inbounds HT .- n_df.HS)
    DST′ = mean(@inbounds n_df.count ./ (n_df.count .- 1) .* DST)
    HT′ = mean(@inbounds n_df.HS .+ DST′)

    return Dict{Symbol, Float64}(
        :FST => DST / HT,
        :FST′ => DST′ / HT′
    )
end


function f_stat_p(data::PopData; nperm::Int = 1000)
    observed = FST_global(data)
    popnames = data.meta.population    
    perm_dict = Dict{Symbol,Vector{Float64}}(
        :FST => Vector{Float64}(undef, nperm-1),
        :FST′ => Vector{Float64}(undef, nperm-1)
        )
    p = Progress(nperm-1, dt = 1, color = :blue)
    @inbounds @sync for n in 1:nperm-1
        Base.Threads.@spawn begin
            tmp = FST_global(permute_samples!(copy(data.loci), popnames))
            perm_dict[:FST][n] = tmp[:FST]
            perm_dict[:FST′][n] = tmp[:FST′]
            next!(p)
        end
    end
    @inbounds for (k,v) in perm_dict
        observed[k] = (1 + sum(observed[k] .<= v))/nperm
    end

    return observed
end
