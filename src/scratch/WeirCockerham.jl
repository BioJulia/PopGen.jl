function _pairwise_wc(data::PopData)
    loc_names = loci(data)
    n_loci = length(loc_names)
    n_pops = length(unique(data.meta.population))
    # get genotype counts
    per_locpop = groupby(data.loci, [:locus, :population])
    nonmissing_counts = DataFrames.combine(per_locpop, :genotype => nonmissing => :n)
    # reshape as a matrix of loci x pop (row x col)
    # collect?
    n_per_locpop = reshape(nonmissing_counts.n', (n_pops, n_loci))'
    n_total = reduce(+, eachcol(n_per_locpop))
    # screen for completely absent loci
    missing_loc = findall(iszero, n_total)
    present_loc = findall(!iszero, n_total)
    if !isempty(missing_loc)
        # index locus names by the missing indices
        tmp_data = exclude(data, locus = loc_names[missing_loc])
        loc_names = loci(tmp_data)
        n_loci = length(loc_names)
        per_locpop = groupby(tmp_data.loci, [:locus, :population])
        nonmissing_counts = DataFrames.combine(per_locpop, :genotype => nonmissing => :n)
        # reshape as a matrix of loci x pop (row x col)
        # collect?
        n_per_locpop = reshape(nonmissing_counts.n', (n_pops, length(present_loci)))'
        n_total = reduce(+, eachcol(n_per_locpop))
    else
        tmp_data = data
    end
    n_pop_per_loc = map(count_nonzeros, eachrow(n_per_locpop))
    per_loc = groupby(tmp_data.loci, :locus)
    allele_counts = DataFrames.combine(
        per_loc, 
        :genotype => allele_count => :allele_count,  # allele counts
        :genotype => allele_freq => :freqs,          # global allele freqs
    )
    # allele freqs per locus per population
    alle_freqs = DataFrames.combine(per_locpop, :genotype => allele_freq => :freqs)
    n_expanded = reduce(hcat, repeat.(eachrow(n_per_locpop), 1, allele_counts.allele_count)) |> permutedims
    corr_n_per_loc = (n_total .- (sum.(eachrow(n_per_locpop .^2)) ./ n_total)) ./ (n_pop_per_loc .- 1) 
    n_tot_expanded = reduce(vcat, repeat.(eachrow(n_total), allele_counts.allele_count))
    corr_n_per_loc_exp = reduce(vcat, repeat.(eachrow(corr_n_per_loc), allele_counts.allele_count))
    _alleles = Base.Iterators.flatten(allele_counts.freqs)  |> collect
    #global_freq = Base.Iterators.flatten(allele_counts.freqs)

end


