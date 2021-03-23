## Nei 1987 FST ##
function nei_fst(data::AbstractDataFrame)
    # observed/expected het per locus per pop
    het_df = DataFrames.combine(
        groupby(data, [:locus, :population]),
        :genotype => nonmissing => :n,
        :genotype => hetero_o => :het_obs,
        :genotype => hetero_e => :het_exp,
        :genotype => allele_freq => :alleles
    )
    # collapse down to retrieve averages and counts
    n_df = DataFrames.combine(
        groupby(het_df, :locus),
        :n => count_nonzeros => :count,
        :n => (n -> count_nonzeros(n) / reciprocal_sum(n)) => :mn,
        [:het_obs, :het_exp, :n] => ((o,e,n) -> mean(skipmissing(gene_diversity_nei87.(e, o, count_nonzeros(n) / reciprocal_sum(n))))) => :HS,
        :het_obs => (o -> mean(skipmissing(o)))=> :Het_obs,
        :alleles => (alleles ->  sum(values(avg_allele_freq(alleles, 2))))=> :avg_freq
        )

    Ht = @inbounds 1.0 .- n_df.avg_freq .+ (n_df.HS ./ n_df.mn ./ n_df.count) - (n_df.Het_obs ./ 2.0 ./ n_df.mn ./ n_df.count)
    DST = @inbounds Ht .- n_df.HS
    DST′ = @inbounds n_df.count ./ (n_df.count .- 1) .* DST
    HT′ = @inbounds n_df.HS .+ DST′
    #TODO replace safemean
    round(mean(skipinfnan(DST′)) / mean(skipinfnan(HT′)), digits = 5)
end

function _pairwise_Nei(data::PopData)
    idx_pdata = groupby(data.loci, :population)
    pops = getindex.(keys(idx_pdata), :population)
    npops = length(idx_pdata)
    results = zeros(npops, npops)
    @sync for i in 2:npops
        Base.Threads.@spawn begin
            for j in 1:(i-1)
                results[i,j] = nei_fst(DataFrame(idx_pdata[[j,i]]))
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
                results[i,j] = weircockerham_fst2(pop1,pop2)
           end
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir & Cockerham")
end