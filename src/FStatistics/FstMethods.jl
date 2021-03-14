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
function _wc_helper(x::AbstractArray{Float64,2})
    view(x,:,1) ./ sum(x, dims = 2)
end

function weircockerham_fst(data::AbstractDataFrame)
    loc_names = unique(data.locus)
    n_loci = length(loc_names)
    n_pops = length(unique(data.population))
    # get genotype counts
    per_locpop = groupby(data, [:locus, :population])
    nonmissing_counts = DataFrames.combine(per_locpop, :genotype => nonmissing => :n)
    # reshape as a matrix of loci x pop (row x col)
    n_per_locpop = reshape(nonmissing_counts.n', (n_pops, n_loci))'
    n_total = reduce(+, eachcol(n_per_locpop))
    # screen for completely absent loci
    missing_loc = findall(iszero, n_total)
    if !isempty(missing_loc)
        present_loc = n_loci - length(missing_loc)
        # index locus names by the missing indices
        tmp_data = filter(:locus => loc -> loc ∉ loc_names[missing_loc], data)
        loc_names = unique(tmp_data.locus)
        n_loci = length(loc_names)
        per_locpop = groupby(tmp_data, [:locus, :population])
        nonmissing_counts = DataFrames.combine(per_locpop, :genotype => nonmissing => :n)
        # reshape as a matrix of loci x pop (row x col)
        n_per_locpop = reshape(nonmissing_counts.n', (n_pops, present_loci))'
        n_total = reduce(+, eachcol(n_per_locpop))
    else
        tmp_data = data
    end

    n_pop_per_loc = map(count_nonzeros, eachrow(n_per_locpop))
    per_loc = groupby(tmp_data, :locus)
    # find all the alleles for
    allele_counts = DataFrames.combine(
        per_loc, 
        :genotype => allele_count => :allele_count,  # allele counts (global)
        :genotype => allele_freq => :freqs,          # allele freqs (global)
    )
    # allele freqs per locus per population
    alle_freqs = DataFrames.combine(per_locpop, :genotype => allele_freq => :freqs)
    # expand out the n matrix to be same dimensions as unique_alleles x pop
    n_expanded = reduce(hcat, repeat.(eachrow(n_per_locpop), 1, allele_counts.allele_count)) |> permutedims
    # expand n_total matrix to be same dimensions as unique_alleles x pop
    n_tot_expanded = reduce(vcat, repeat.(eachrow(n_total), allele_counts.allele_count))
    # calculate corrected n per locus
    corr_n_per_loc = (n_total .- (sum(n_per_locpop .^2, dims = 2) ./ n_total)) ./ (n_pop_per_loc .- 1) 
    # expand corr_n matrix to be same dimensions as unique_alleles x pop
    corr_n_per_loc_exp = reduce(vcat, repeat.(eachrow(corr_n_per_loc), allele_counts.allele_count))
    # list of alleles in each locus
    _alleles_perloc = [sort(unique_alleles(i.genotype)) for i in per_loc]
    # extremely convoluted, creates a big matrix of allele freqs per locus per population
    # TODO are there too many reshapes going on?
    #TODO move into its own function? This seems like it could be a recurring thing
    afreq_tmp = permutedims(reshape(alle_freqs.freqs, n_pops, :))
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
    _freqs = map(i -> values(sort(i)), allele_counts.freqs) |> Base.Iterators.flatten |> collect

    #heterozygotes per allele per locus per population
    # gets emptied from the popfirst! calls below(!)
    genos_vec = [i.genotype for i in per_locpop]

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
    n_correction = reduce(vcat, fill.(n_pop_per_loc, allele_counts.allele_count))

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
    results = zeros(npops, npops)
    @sync for i in 2:npops
        Base.Threads.@spawn begin
            for j in 1:(i-1)
                results[i,j] = weircockerham_fst(DataFrame(idx_pdata[[j,i]]))
            end
        end
    end
    return PairwiseFST(DataFrame(results, Symbol.(pops)), "Weir & Cockerham")
end