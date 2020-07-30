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


function bootstrap_relatedness_before(data::PopData, ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; method::F, iterations::Int) where T <: GenoArray where U <: NamedTuple where F
    loci_gdf = groupby(data.loci, :locus)
    relate_vec_boot = Vector{Union{Missing,Float64}}(undef, iterations)
    n_loc = length(locus_names)
    for iter in 1:iterations
        # bootstrap the indices
        boot_idx = rand(1:n_loc, n_loc)
        # sample the source vectors with the resampled/bootstrapped indices
        ind1_boot, ind2_boot, loc_boot = map(i -> getindex(i, boot_idx), [ind1, ind2, locus_names])
        #n_per_loci = map(i -> nonmissing(data, i), loc_boot)
        # faster/cheaper n counting
        n_per_loci = map(i -> nonmissing(loci_gdf[(i,)].genotype, loc_boot))
        loc_samp,gen_samp1,gen_samp2, n_per_loc = collect.(skipmissings(Symbol.(loc_boot), ind1_boot, ind2_boot, n_per_loci))

        relate_vec_boot[iter] = method(gen_samp1, gen_samp2, loc_samp, alleles, loc_n = n_per_loc, n_samples = n_loc)
    end
    return relate_vec_boot
end

## bootstrap after removing missing
function bootstrap_relatedness_after(data::PopData, ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; method::F, iterations::Int) where T <: GenoArray where U <: NamedTuple where F
    loci_gdf = groupby(data.loci, :locus)
    relate_vec_boot = Vector{Union{Missing,Float64}}(undef, iterations)
    shared_loci, ind1_geno, ind2_geno = collect.(skipmissings(Symbol.(loc_boot), ind1_boot, ind2_boot))
    n_loci = length(shared_loci)
    for iter in 1:iterations
        # bootstrap the indices
        boot_idx = rand(1:n_loc, n_loc)
        # sample the source vectors with the resampled/bootstrapped indices
        ind1_boot, ind2_boot, loc_boot = map(i -> getindex(i, boot_idx), [ind1_geno, ind2_geno, shared_loci])
        #n_per_loci = map(i -> nonmissing(data, i), loc_boot)
        # faster/cheaper n counting
        n_per_loci = map(i -> nonmissing(loci_gdf[(i,)].genotype, loc_boot))

        relate_vec_boot[iter] = method(gen_samp1, gen_samp2, loc_samp, alleles, loc_n = n_per_loc, n_samples = n_loci)
    end
    return relate_vec_boot
end