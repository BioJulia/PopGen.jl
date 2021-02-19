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
    _alleles = map(keys, allele_counts.freqs) |> Base.Iterators.flatten |> collect
    _freqs = map(values, allele_counts.freqs) |> Base.Iterators.flatten |> collect

    # # heterozygotes per allele per locus per population
    #TODO remove sorting (it's there for testing)
    # sorted list of alleles in each locus
    _alleles_perloc = map(j -> sort(collect(keys(j))), allele_counts.freqs)
    
    genos_vec = [i.genotype for i in per_locpop]

    # create matrix of heterozygote occurrences per allele per pop
    het_mtx = reduce(vcat,
        # map across the vector of alleles for each locus
        # vcat will concatenate the returned matrices into a single matrix
        map(_alleles_perloc) do _alleles
            reduce(hcat,
                # each element in x is locus × population, so use a comprehension to
                # do counthet() as many times as there are populations, popfirst!'ing
                # the first element of x until it's ultimately empty
                # then concatenate it into a matrix
                [counthet(popfirst!(genos_vec), _alleles) for pop_n in 1: n_pops]
            )
        end
    )
end

x = [i.genotype for i in zz] ;
ns = length.(_alleles_perloc)


tst = reduce(vcat,
        # map across the vector of alleles for each locus
        # vcat will concatenate the returned matrices into a single matrix
        map(_alleles_perloc) do _alleles
            reduce(hcat,
                # each element in x is locus × population, so use a comprehension to
                # do counthet() as many times as there are populations, popfirst!'ing
                # the first element of x until it's ultimately empty
                # then concatenate it into a matrix
                [counthet(popfirst!(x), _alleles) for pop_n in 1:pops]
            )
        end
    )