"""
    counthet(geno::T, allele::Int) where T<:GenoArray
Given a `GenoArray`, count the number of times `allele` appears in the
heterozygous state.
"""
function counthet(geno::T, allele::Int) where T<:GenoArray
    reduce(+, (allele .∈ skipmissing(geno)) .& (ishet(skipmissing(geno))))
end


"""
    counthom(geno::T, allele::Int) where T<:GenoArray
Given a `GenoArray`, count the number of times `allele` appears in the
homozygous state.
"""
function counthom(geno::T, allele::Int) where T<:GenoArray
    reduce(+, (allele .∈ skipmissing(geno)) .& (ishom(skipmissing(geno))))
end


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
    _alleles_perloc = map(keys, allele_counts.freqs)
    for i in groupby(tmp_data.loci, [:locus, :population])
        ishet(i.genotype)
    end
    #return groupby(tmp_data.loci, [:locus, :population])
end


