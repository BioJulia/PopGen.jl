function _pairwise_Hudson_lxl(data::PopData)
    !isbiallelic(data) && throw(error("Data must be biallelic to use the Hudson estimator"))
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    nloci = data.metadata.loci
    locnames = loci(data)
    allpairs = pairwisepairs(pops)
    npairs = length(collect(allpairs))
    results = Vector{Vector{Float64}}(undef, npairs)
    p1 = repeat(getindex.(allpairs,1), inner = nloci)
    p2 = repeat(getindex.(allpairs,2), inner = nloci)
    locs = repeat(locnames, outer = npairs)
    @inbounds @sync for (i,j) in enumerate(allpairs)
        Base.Threads.@spawn begin
            pop1 = reshape(idx_pdata[(population = j[1],)].genotype, :, nloci)
            pop2 = reshape(idx_pdata[(population = j[2],)].genotype, :, nloci)
            results[i] = _hudson_fst_lxl(pop1, pop2)
        end
    end
    return DataFrame(:pop1 => PooledArray(p1, compress = true), :pop2 => PooledArray(p2, compress = true),:locus => PooledArray(locs, compress = true) ,:fst => reduce(vcat, results))
end

function _hudson_fst_lxl(population_1::T, population_2::T) where T<:AbstractMatrix
    @views [_hudson_fst(population_1[:,i], population_2[:,i]) for i in 1:size(population_1,2)]
end

function _nei_fst_lxl(population_1::T, population_2::T) where T<:AbstractMatrix
    n =  hcat(map(nonmissing, eachcol(population_1)), map(nonmissing, eachcol(population_2)))
    # number of populations represented per locus
    n_pop_per_loc = map(countnonzeros, eachrow(n)) 
    # corrected n for population size
    corr_n_per_loc = n_pop_per_loc ./ map(reciprocalsum, eachrow(n))
    # observed heterozygosity
    het_obs_p1 = map(_hetero_obs, eachcol(population_1))
    het_obs_p2 = map(_hetero_obs, eachcol(population_2))
    # expected heterozygosity
    het_exp_p1 = map(_hetero_exp, eachcol(population_1))
    het_exp_p2 = map(_hetero_exp, eachcol(population_2))
    # genic diversity for population 1 and 2
    p1_nei = _genediversitynei87.(het_exp_p1, het_obs_p1, corr_n_per_loc)
    p2_nei = _genediversitynei87.(het_exp_p2, het_obs_p2, corr_n_per_loc)
    # mean genic diversity
    HS = map(i -> mean(skipmissing(i)), zip(p1_nei, p2_nei))
    # allele freqs for population 1 and 2
    alle_freq_p1 = map(allelefreq, eachcol(population_1))
    alle_freq_p2 = map(allelefreq, eachcol(population_2))
    # sum of the squared average allele frequencies
    avg_freq = avg_allelefreq.(zip(alle_freq_p1, alle_freq_p2),2) .|> values .|> sum
    # average observed heterozygosity per locus
    Het_obs = map(i -> mean(skipmissing(i)), zip(het_obs_p1, het_obs_p2))

    Ht = @inbounds 1.0 .- avg_freq .+ (HS ./ corr_n_per_loc ./ n_pop_per_loc) - (Het_obs ./ 2.0 ./ corr_n_per_loc ./ n_pop_per_loc)
    DST = @inbounds Ht .- HS
    DST′ = @inbounds n_pop_per_loc ./ (n_pop_per_loc .- 1) .* DST
    HT′ = @inbounds HS .+ DST′
    round.(skipinfnan(DST′) ./ skipinfnan(HT′), digits = 5)
end

function _pairwise_Nei_lxl(data::PopData)
    idx_pdata = groupby(data.genodata, :population)
    pops = getindex.(keys(idx_pdata), :population)
    nloci = data.metadata.loci
    locnames = loci(data)
    allpairs = pairwisepairs(pops)
    npairs = length(collect(allpairs))
    results = Vector{Vector{Float64}}(undef, npairs)
    p1 = repeat(getindex.(allpairs,1), inner = nloci)
    p2 = repeat(getindex.(allpairs,2), inner = nloci)
    locs = repeat(locnames, outer = npairs)
    @inbounds @sync for (i,j) in enumerate(allpairs)
        Base.Threads.@spawn begin
            pop1 = reshape(idx_pdata[(population = j[1],)].genotype, :, nloci)
            pop2 = reshape(idx_pdata[(population = j[2],)].genotype, :, nloci)
            results[i] = _nei_fst_lxl(pop1, pop2)
        end
    end
    return DataFrame(:pop1 => PooledArray(p1, compress = true), :pop2 => PooledArray(p2, compress = true),:locus => PooledArray(locs, compress = true) ,:fst => reduce(vcat, results))
end