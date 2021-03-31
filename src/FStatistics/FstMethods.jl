function _pairwise_Hudson(data::PopData)
    !isbiallelic(data) && throw(error("Data must be biallelic to use the Hudson estimator"))
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    n_loci = size(data)[2]
    results = zeros(npops, npops)
    @sync for i in 2:npops
        Base.Threads.@spawn begin
            for j in 1:(i-1)
                pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
                pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
                results[i,j] = fst_hudson(pop1,pop2)
           end
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Hudson et al. 1992")
end

function hudson_fst(population_1::T, population_2::T) where T<:AbstractMatrix
    fst_perloc = @views [_hudson_fst(population_1[:,i], population_2[:,i]) for i in 1:size(population_1,2)]
    return mean(skipmissing(skipinfnan(fst_perloc)))
end

# helper function to do the math for Hudson FST on a locus
function _hudson_fst(pop1::T, pop2::T) where T<:GenoArray
    p1_frq = allele_freq(pop1)
    p2_frq = allele_freq(pop2)
    # find the shared allele(s) and choose one of them to be "P"
    # this is a safeguard if one population is completely homozygous for an allele
    p_allele = intersect(keys(p1_frq), keys(p2_frq)) |> first
    p1 = p1_frq[p_allele]
    q1 = 1.0 - p1
    p2 = p2_frq[p_allele]
    q2 = 1.0 - p2 
    # degrees of freedom is the number of alleles - 1
    df1 = (count(!ismissing, pop1) * 2) - 1
    df2 = (count(!ismissing, pop2) * 2) - 1
    numerator = (p1 - p2)^2 - (p1*q1/df1) - (p2*q2/df2)
    denominator = (p1*q2) + (p2*q1)
    return numerator/denominator
end



## Nei 1987 FST ##
function nei_fst(population_1::T, population_2::T) where T<:AbstractMatrix
    n =  hcat(map(nonmissing, eachcol(population_1)), map(nonmissing, eachcol(population_2)))
    # number of populations represented per locus
    n_pop_per_loc = map(count_nonzeros, eachrow(n)) 
    # corrected n for population size
    corr_n_per_loc = n_pop_per_loc ./ map(reciprocal_sum, eachrow(n))
    # observed heterozygosity
    het_obs_p1 = map(hetero_o, eachcol(population_1))
    het_obs_p2 = map(hetero_o, eachcol(population_2))
    # expected heterozygosity
    het_exp_p1 = map(hetero_e, eachcol(population_1))
    het_exp_p2 = map(hetero_e, eachcol(population_2))
    # genic diversity for population 1 and 2
    p1_nei = gene_diversity_nei87.(het_exp_p1, het_obs_p1, corr_n_per_loc)
    p2_nei = gene_diversity_nei87.(het_exp_p2, het_obs_p2, corr_n_per_loc)
    # mean genic diversity
    HS = map(i -> mean(skipmissing(i)), zip(p1_nei, p2_nei))
    # allele freqs for population 1 and 2
    alle_freq_p1 = map(allele_freq, eachcol(population_1))
    alle_freq_p2 = map(allele_freq, eachcol(population_2))
    # sum of the squared average allele frequencies
    avg_freq = avg_allele_freq.(zip(alle_freq_p1, alle_freq_p2),2) .|> values .|> sum
    # average observed heterozygosity per locus
    Het_obs = map(i -> mean(skipmissing(i)), zip(het_obs_p1, het_obs_p2))

    Ht = @inbounds 1.0 .- avg_freq .+ (HS ./ corr_n_per_loc ./ n_pop_per_loc) - (Het_obs ./ 2.0 ./ corr_n_per_loc ./ n_pop_per_loc)
    DST = @inbounds Ht .- HS
    DST′ = @inbounds n_pop_per_loc ./ (n_pop_per_loc .- 1) .* DST
    HT′ = @inbounds HS .+ DST′
    round(mean(skipinfnan(DST′)) / mean(skipinfnan(HT′)), digits = 5)
end

function _pairwise_Nei(data::PopData)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    n_loci = size(data)[2]
    results = zeros(npops, npops)
    @sync for i in 2:npops
        Base.Threads.@spawn begin
            for j in 1:(i-1)
                pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
                pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
                results[i,j] = nei_fst(pop1,pop2)
           end
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Nei 1987")
end


## Weir & Cockerham 1984 FST ##
function weircockerham_fst(population_1::T, population_2::T) where T<:AbstractMatrix
    n_loci = size(population_1, 2)
    n_pops = 2
    # get genotype counts
    # reshape as a matrix of loci x pop (row x col)
    n_per_locpop =  hcat(map(nonmissing, eachcol(population_1)), map(nonmissing, eachcol(population_2)))
    n_total = sum(n_per_locpop, dims = 2)
    # screen for completely absent loci
    present_loc = 0 .∉ eachrow(n_per_locpop)
    if 0 ∈ present_loc
        # index locus names by the missing indices
        pop_1 = @view population_1[:,present_loc]
        pop_2 = @view population_2[:,present_loc]
        n_per_locpop =  hcat(map(nonmissing, eachcol(pop_1)), map(nonmissing, eachcol(pop_2)))
        n_total = sum(n_per_locpop, dims = 2)
    else
        pop_1 = population_1
        pop_2 = population_2
    end
    merged = vcat(pop_1, pop_2)
    n_pop_per_loc = map(count_nonzeros, eachrow(n_per_locpop))

    # global allele counts
    glob_allele_counts = map(allele_count, eachcol(merged))
    # global allele frequencies
    glob_allele_freqs = map(allele_freq, eachcol(merged))
    
    # allele freqs per locus per population
    pop_1_freq = map(allele_freq, eachcol(pop_1))
    pop_2_freq = map(allele_freq, eachcol(pop_2))
    #alle_freqs = DataFrames.combine(per_locpop, :genotype => allele_freq => :freqs)
    # expand out the n matrix to be same dimensions as unique_alleles x pop
    n_expanded = reduce(hcat, repeat.(eachrow(n_per_locpop), 1, glob_allele_counts)) |> permutedims
    # expand n_total matrix to be same dimensions as unique_alleles x pop
    n_tot_expanded = reduce(vcat, repeat.(eachrow(n_total), glob_allele_counts))
    # calculate corrected n per locus
    corr_n_per_loc = (n_total .- (sum(n_per_locpop .^2, dims = 2) ./ n_total)) ./ (n_pop_per_loc .- 1) 
    # expand corr_n matrix to be same dimensions as unique_alleles x pop
    corr_n_per_loc_exp = reduce(vcat, repeat.(eachrow(corr_n_per_loc), glob_allele_counts))
    # list of alleles in each locus
    _alleles_perloc = [sort(unique_alleles(i)) for i in eachcol(merged)]
    # extremely convoluted, creates a big matrix of allele freqs per locus per population
    # TODO are there too many reshapes going on?
    #TODO move into its own function? This seems like it could be a recurring thing
    afreq_tmp = hcat(pop_1_freq, pop_2_freq)
    allele_freq_pop = reshape(
        reduce(vcat,
            map(zip(eachrow(afreq_tmp), _alleles_perloc)) do (_freqs_p, _alle)
                reduce(hcat,
                    map(_freqs_p) do _freqs
                        [get(_freqs, i, 0.) for i in _alle]     # query the dict of allele freqs
                    end
                )
            end
            ),
       :, n_pops  # reshape by these dimensions
    )
    # global allele freqs
    _freqs = map(i -> values(sort(i)), glob_allele_freqs) |> Base.Iterators.flatten |> collect

    #heterozygotes per allele per locus per population
    # gets emptied from the popfirst! calls below(!)
    _p1 = collect(eachcol(pop_1))
    _p2 = collect(eachcol(pop_2))
    genos_vec = permutedims([_p1 _p2])[:]
    # create matrix of heterozygote occurrences per allele per pop
    het_mtx = reduce(vcat,     # vcat will concatenate the returned matrices into a single matrix
        map(_alleles_perloc) do _alleles          # map across the vector of alleles for each locus
            reduce(hcat,
                # each element in x is locus × population, so use a comprehension to
                # do counthet() as many times as there are populations, popfirst!'ing
                # the first element of x until it's ultimately empty
                # then concatenate it into a matrix
                [counthet(popfirst!(genos_vec), _alleles) for pop_n in 1: n_pops]
            )
        end
    )
    μ_het =  (2 * n_expanded .* allele_freq_pop - het_mtx) / 2
    SSG = sum(n_expanded .* allele_freq_pop - μ_het, dims = 2)
    SSi = sum(n_expanded .* (allele_freq_pop - 2 * allele_freq_pop .^ 2) + μ_het, dims = 2)
    SSP = 2 .* sum(n_expanded .* reduce(hcat, map(i -> (i .- _freqs) .^ 2, eachcol(allele_freq_pop))), dims = 2)
    n_correction = reduce(vcat, fill.(n_pop_per_loc, glob_allele_counts))

    MSG = SSG ./ n_tot_expanded
    MSP = SSP ./ (n_correction .- 1)
    MSI = SSi ./ (n_tot_expanded - n_correction)
    σ_w = MSG
    σ_b = 0.5 * (MSI - MSG)
    σ_a = 0.5 ./ corr_n_per_loc_exp .* (MSP - MSI)
    σ = hcat(σ_a, σ_b, σ_w)
    σ_total = map(i -> sum(skipinfnan(i)), eachcol(σ))
    fst_total = round(σ_total[1] / sum(σ_total), digits = 5)
end


function _pairwise_WeirCockerham(data::PopData)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    n_loci = size(data)[2]
    results = zeros(npops, npops)
    @sync for i in 2:npops
        Base.Threads.@spawn begin
            for j in 1:(i-1)
                pop1 = reshape(idx_pdata[i].genotype, :, n_loci)
                pop2 = reshape(idx_pdata[j].genotype, :, n_loci)
                results[i,j] = weircockerham_fst(pop1,pop2)
           end
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir & Cockerham")
end
