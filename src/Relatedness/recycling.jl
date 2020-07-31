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
        loc_samp,gen_samp1,gen_samp2, n_per_loc = collect.(skipmissings(loc_boot, ind1_boot, ind2_boot, n_per_loci))

        relate_vec_boot[iter] = method(gen_samp1, gen_samp2, loc_samp, alleles, loc_n = n_per_loc, n_samples = n_loc)
    end
    return relate_vec_boot
end

## bootstrap after removing missing
## replace data::PopData with data::GroupDataFrame?
function bootstrap_relatedness_after(data::PopData, ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; method::F, iterations::Int) where T <: GenoArray where U <: NamedTuple where F
    loci_gdf = groupby(data.loci, :locus)
    relate_vec_boot = Vector{Union{Missing,Float64}}(undef, iterations)
    n_loc = length(locus_names)
    for iter in 1:iterations
        # bootstrap the indices
        boot_idx = rand(1:n_loc, n_loc)
        # sample the source vectors with the resampled/bootstrapped indices
        ind1_boot, ind2_boot, loc_boot = map(i -> getindex(i, boot_idx), [ind1, ind2, locus_names])
        # faster/cheaper n counting
        n_per_loci = map(i -> nonmissing(loci_gdf[(i,)].genotype, loc_boot))

        relate_vec_boot[iter] = method(ind1_boot, ind2_boot, loc_boot, alleles, loc_n = n_per_loci, n_samples = n_loc)
    end
    return relate_vec_boot
end



function relatedness_bootstrap_before(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975)) where F
    loci_names = Symbol.(loci(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_samples = length(samples(data))
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    boot_means, boot_medians, boot_ses = map(i -> deepcopy(relate_vecs), 1:3)
    boot_CI = map(i -> Vector{Union{Missing,Tuple{Float64,Float64}}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1

            geno2 = popdata_idx[(ind2,)].genotype

            # filter out loci missing in at least one individual in the pair
            loc, gen1, gen2, n_per_loc = collect.(skipmissings(loci_names, geno1, geno2, n_per_loci))

            # populate shared_loci array
            @inbounds shared_loci[idx] = length(loc)

            @inbounds for (i, mthd) in enumerate(method)
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples)
                boot_out = bootstrap_locus_before(data, mthd, geno1, geno2, loci_names, allele_frequencies, method = mthd, iterations = iterations)
                @inbounds boot_means[i][idx], boot_medians[i][idx], boot_ses[i][idx], boot_CI[i][idx] = bootstrap_summary(boot_out, iterations, interval)
            end
            update!(p, idx)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    boot_mean_colnames = [Symbol("$i"*"_mean") for i in method]
    boot_median_colnames = [Symbol("$i"*"_median") for i in method]
    boot_se_colnames = [Symbol("$i"*"_SE") for i in method]
    CI_percent = convert(Int64, round(interval[2] - interval[1], digits = 2) * 100)
    boot_CI_colnames = [Symbol("$i"*"_CI_"*"$CI_percent") for i in method]

    out_df = DataFrame(:sample_1 => map(i -> i[1], sample_pairs), :sample_2 => map(i -> i[2], sample_pairs), :n_loci => shared_loci)
    @inbounds for (i, mth) in enumerate(method_colnames)
        out_df[:, mth] = relate_vecs[i]
        out_df[:, boot_mean_colnames[i]] = boot_means[i]
        out_df[:, boot_median_colnames[i]] = boot_medians[i]
        out_df[:, boot_se_colnames[i]] = boot_ses[i]
        out_df[:, boot_CI_colnames[i]] = boot_CI[i]
    end

    return out_df
end

function relatedness_bootstrap_after(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975)) where F
    loci_names = Symbol.(loci(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_samples = length(samples(data))
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    boot_means, boot_medians, boot_ses = map(i -> deepcopy(relate_vecs), 1:3)
    boot_CI = map(i -> Vector{Union{Missing,Tuple{Float64,Float64}}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1

            geno2 = popdata_idx[(ind2,)].genotype

            # filter out loci missing in at least one individual in the pair
            loc, gen1, gen2, n_per_loc = collect.(skipmissings(loci_names, geno1, geno2, n_per_loci))

            # populate shared_loci array
            @inbounds shared_loci[idx] = length(loc)
            
            @inbounds for (i, mthd) in enumerate(method)
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples)
                boot_out = bootstrap_locus_before(data, mthd, gen1, gen2, loc, allele_frequencies, method = mthd, iterations = iterations)
                @inbounds boot_means[i][idx], boot_medians[i][idx], boot_ses[i][idx], boot_CI[i][idx] = bootstrap_summary(boot_out, iterations, interval)
            end
            update!(p, idx)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    boot_mean_colnames = [Symbol("$i"*"_mean") for i in method]
    boot_median_colnames = [Symbol("$i"*"_median") for i in method]
    boot_se_colnames = [Symbol("$i"*"_SE") for i in method]
    CI_percent = convert(Int64, round(interval[2] - interval[1], digits = 2) * 100)
    boot_CI_colnames = [Symbol("$i"*"_CI_"*"$CI_percent") for i in method]

    out_df = DataFrame(:sample_1 => map(i -> i[1], sample_pairs), :sample_2 => map(i -> i[2], sample_pairs), :n_loci => shared_loci)
    @inbounds for (i, mth) in enumerate(method_colnames)
        out_df[:, mth] = relate_vecs[i]
        out_df[:, boot_mean_colnames[i]] = boot_means[i]
        out_df[:, boot_median_colnames[i]] = boot_medians[i]
        out_df[:, boot_se_colnames[i]] = boot_ses[i]
        out_df[:, boot_CI_colnames[i]] = boot_CI[i]
    end

    return out_df
end


















# extremely incomplete
function relatedness_bootstrap(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975)) where F
    loci_names = Symbol.(loci(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_loci = length(loci_names)
    n_samples = length(samples(data))
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    boot_means, boot_medians, boot_ses = map(i -> deepcopy(relate_vecs), 1:3)
    boot_CI = map(i -> Vector{Union{Missing,Tuple{Float64,Float64}}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    for i in 1:iterations
        boot_idx = rand(1:n_loci, n_loci)
        @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
            geno1 = popdata_idx[(ind1,)].genotype
            @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
                idx += 1
                geno2 = popdata_idx[(ind2,)].genotype

                ind1_boot, ind2_boot, loc_boot, n_per_loc = map(i -> getindex(i, boot_idx), [ind1, ind2, loci_names, n_per_loc])
                # faster/cheaper n counting
                #loc, gen1, gen2, n_per_loc_sm = collect.(skipmissings(loc_boot, ind1_boot, ind2_boot, n_per_loci))
                # filter out loci missing in at least one individual in the pair
                loc, gen1, gen2, n_per_loc_sm = collect.(skipmissings(loci_names[boot_idx], geno1[boot_idx], geno2[boot_idx], n_per_loci[boot_idx]))

                # populate shared_loci array
                @inbounds shared_loci[idx] = length(loc)
                @inbounds for (i, mthd) in enumerate(method)
                    @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples)
                    boot_out = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples)
                    #boot_out = bootstrap_locus(data, mthd, ind1, ind2, iterations, allele_frequencies)
                    @inbounds boot_means[i][idx], boot_medians[i][idx], boot_ses[i][idx], boot_CI[i][idx] = bootstrap_summary(boot_out, iterations, interval)
                end
            update!(p, idx)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    boot_mean_colnames = [Symbol("$i"*"_mean") for i in method]
    boot_median_colnames = [Symbol("$i"*"_median") for i in method]
    boot_se_colnames = [Symbol("$i"*"_SE") for i in method]
    CI_percent = convert(Int64, round(interval[2] - interval[1], digits = 2) * 100)
    boot_CI_colnames = [Symbol("$i"*"_CI_"*"$CI_percent") for i in method]

    out_df = DataFrame(:sample_1 => map(i -> i[1], sample_pairs), :sample_2 => map(i -> i[2], sample_pairs), :n_loci => shared_loci)
    @inbounds for (i, mth) in enumerate(method_colnames)
        out_df[:, mth] = relate_vecs[i]
        out_df[:, boot_mean_colnames[i]] = boot_means[i]
        out_df[:, boot_median_colnames[i]] = boot_medians[i]
        out_df[:, boot_se_colnames[i]] = boot_ses[i]
        out_df[:, boot_CI_colnames[i]] = boot_CI[i]
    end

    return out_df
end