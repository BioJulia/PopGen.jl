# matrix version?
function QuellerGoodnight2(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple

    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    numerator2 = 0.0
    denominator1 = 0.0
    denominator2 = 0.0

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    # drop missing values
    loc, geno1, geno2 = collect.(skipmissings(Symbol.(loci(data)), geno1, geno2))

    # convert to 2-dimensional matrix for each individual
    gen1 = hcat(getindex.(geno1, 1), getindex.(geno1, 2))
    gen2 = hcat(getindex.(geno2, 1), getindex.(geno2, 2))
    # ((a == c) + (a == d) + (b == c) + (b == d))
    id_mtx = reduce(*, (gen1 .== gen2) + (gen1 .== reverse(gen2, dims = 2)), dims = 2)
    a_b_frqs = [alleles[loc[i]][gen1[i,j]] for i in 1:size(gen1)[1], j in 1:2]
    c_d_frqs = [alleles[loc[i]][gen2[i,j]] for i in 1:size(gen2)[1], j in 1:2]
    
    numerator1 = reduce(+, id_mtx - (2 * reduce(+, a_b_frqs, dims = 2))) 
    numerator2 = reduce(+, id_mtx - (2 * reduce(+, c_d_frqs, dims = 2)))
    # reduce(+,)
    denominator1 = 2.0 * (1.0 + map(x ->x[1] == x[2],eachrow(gen1)) - foldl(-, eachcol(a_b_frqs)))
    denominator2 = 2.0 * (1.0 + map(x ->x[1] == x[2],eachrow(gen2)) - foldl(-, eachcol(a_b_frqs)))
    #return [alleles[loc[i]][gen1[i,j]] for i in 1:size(gen1)[1], j in 1:2]
    return numerator2
    a,b = gen1
    c,d = gen2

        ident = ((a == c) + (a == d) + (b == c) + (b == d))
        numerator1 += ident - 2.0 * (alleles[loc][a] + alleles[loc][b])
        numerator2 += ident - 2.0 * (alleles[loc][c] + alleles[loc][d])

        denominator1 += (2.0 * (1.0 + (a==b) - alleles[loc][a] - alleles[loc][b]))
        denominator2 += (2.0 * (1.0 + (c==d) - alleles[loc][c] - alleles[loc][d]))

    return (numerator1/denominator1 + numerator2/denominator2)/2.0
end


function relatedness_no_boot_gdf(data::PopData, sample_names::Vector{String}; method::F) where F
    loci_names = Symbol.(loci(data))
    n_samples = length(samples(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data, :name)
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1

            geno2 = popdata_idx[(ind2,)].genotype

            # filter out loci missing in at least one individual in the pair
            loc,gen1,gen2, n_per_loc = collect.(skipmissings(loci_names, geno1, geno2, n_per_loci))

            # populate shared_loci array
            @inbounds shared_loci[idx] = length(loc)
            @inbounds [relate_vecs[i][idx] = mth(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples) for (i,mth) in enumerate(method)]

            update!(p, idx)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci)
    [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    return out_df
end