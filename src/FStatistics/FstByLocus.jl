function _fst_lxl(data::PopData,::Val{:Hudson})
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    nloci = data.metadata.loci
    locnames = loci(data)
    allpairs = pairwisepairs(pops)
    npops = data.metadata.populations
    npairs = Int64(npops * (npops-1) / 2)
    results = Vector{Vector{Union{Missing, Float64}}}(undef, npairs)
    p1 = PooledArray(repeat(getindex.(allpairs,1), inner = nloci), compress = true)
    p2 = PooledArray(repeat(getindex.(allpairs,2), inner = nloci), compress = true)
    locs = PooledArray(repeat(locnames, outer = npairs), compress = true)
    @inbounds @sync for (i,j) in enumerate(allpairs)
        Base.Threads.@spawn begin
            @inbounds pop1 = reshape(idx_pdata[(population = j[1],)].genotype, :, nloci)
            @inbounds pop2 = reshape(idx_pdata[(population = j[2],)].genotype, :, nloci)
            @inbounds results[i] = _hudson_fst_lxl(pop1, pop2)
        end
    end
    return PairwiseFST(
        DataFrame(:population1 => p1, :population2 => p2,:locus => locs, :fst => reduce(vcat, results)),
        "Hudson estimator"
    )
end

function _hudson_fst_lxl(population_1::T, population_2::T) where T<:AbstractMatrix
    @views Union{Float64,Missing}[Hudson(population_1[:,i], population_2[:,i]) for i in 1:size(population_1,2)]
end

function _nei_fst_lxl(population_1::T, population_2::T) where T<:AbstractMatrix
    result = Vector{Union{Missing,Float64}}(undef, size(population_1, 2))
    for (i,p1) in enumerate(eachcol(population_1))
        p2 = view(population_2, :, i)
        n1 = nonmissing(p1)
        n2 = nonmissing(p2)
        # number of populations represented per locus
        n_pop_per_loc = Float64((n1 > 0) + (n2 > 0))
        n_pop_per_loc < 2.0 && continue
        # corrected n for population size
        corr_n_per_loc = n_pop_per_loc / (reciprocal(n1) + reciprocal(n2))
        # observed heterozygosity
        het_obs_p1 = _hetero_obs(p1)
        het_obs_p2 = _hetero_obs(p2)
        # expected heterozygosity
        het_exp_p1 = _hetero_exp(p1)
        het_exp_p2 = _hetero_exp(p2)
        # genic diversity for population 1 and 2
        p1_nei = _genediversitynei87(het_exp_p1, het_obs_p1, corr_n_per_loc)
        p2_nei = _genediversitynei87(het_exp_p2, het_obs_p2, corr_n_per_loc)
        # mean genic diversity
        HS = (p1_nei + p2_nei) / n_pop_per_loc
        alle_freq_p1 = allelefreq(p1)
        alle_freq_p2 = allelefreq(p2)
        avg_freq = sum(values(avg_allelefreq([alle_freq_p1, alle_freq_p2],2)))
        Het_obs = (het_obs_p1 + het_obs_p2) / n_pop_per_loc
        Ht = 1.0 - avg_freq + (HS / corr_n_per_loc / n_pop_per_loc) - (Het_obs / 2.0 / corr_n_per_loc / n_pop_per_loc)
        DST = 2 * (Ht - HS)
        DST′ = DST * (n_pop_per_loc - 1)
        HT′ = HS + DST′
        @inbounds result[i] = round(DST′ / HT′, digits = 5)
  end
  return result
end

function _fst_lxl(data::PopData, ::Val{:Nei})
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(pops)
    nloci = data.metadata.loci
    locnames = loci(data)
    allpairs = pairwisepairs(pops)
    npairs = Int64(npops * (npops-1) / 2)
    results = Vector{Vector{Union{Missing,Float64}}}(undef, npairs)
    @inbounds p1 = PooledArray(repeat(getindex.(allpairs,1), inner = nloci), compress = true)
    @inbounds p2 = PooledArray(repeat(getindex.(allpairs,2), inner = nloci), compress = true)
    locs = PooledArray(repeat(locnames, outer = npairs), compress = true)
    @inbounds @sync for (i,j) in enumerate(allpairs)
        Base.Threads.@spawn begin
            @inbounds pop1 = reshape(idx_pdata[(population = j[1],)].genotype, :, nloci)
            @inbounds pop2 = reshape(idx_pdata[(population = j[2],)].genotype, :, nloci)
            @inbounds results[i] = _nei_fst_lxl(pop1, pop2)
        end
    end
    return PairwiseFST(
        DataFrame(:population1 => p1, :population2 => p2,:locus => locs, :fst => reduce(vcat, results)),
        "Nei estimator"
    )
end