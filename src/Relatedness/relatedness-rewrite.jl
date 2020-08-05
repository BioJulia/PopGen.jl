function relatedness_no_boot_fix2(data::PopData, sample_names::Vector{String}; method::F, inbreeding::Bool) where F
    loci_names = Symbol.(loci(data))
    n_samples = length(samples(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    @inbounds Base.Threads.@threads for i in 1:length(sample_pairs)
        @inbounds sample_1 = sample_pairs[i][1]
        @inbounds sample_2 = sample_pairs[i][2]
        @inbounds geno1 = popdata_idx[(sample_1,)].genotype
        @inbounds geno2 = popdata_idx[(sample_2,)].genotype
        keep_idx = nonmissings(geno1, geno2)
        @inbounds shared_loci[i] = length(keep_idx)
        @inbounds [relate_vecs[j][i] = @inbounds mth(geno1[keep_idx], geno2[keep_idx], loci_names[keep_idx], allele_frequencies, loc_n = n_per_loci[keep_idx], n_samples = n_samples, inbreeding = inbreeding) for (j,mth) in enumerate(method)]
        next!(p)
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci)
    [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    return out_df
end



function relatedness_boot_nonmissing(data::PopData, sample_names::Vector{String}; method::F, iterations::Int, interval::Tuple{Float64, Float64} = (0.025, 0.975), inbreeding::Bool) where F
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
    @inbounds Base.Threads.@threads for i in 1:length(sample_pairs)
        @inbounds sample_1 = sample_pairs[i][1]
        @inbounds sample_2 = sample_pairs[i][2]
        @inbounds geno1 = popdata_idx[(sample_1,)].genotype
        @inbounds geno2 = popdata_idx[(sample_2,)].genotype
        # get index for genotype appearing missing in at least one individual in the pair
        keep_idx = nonmissings(geno1, geno2)
        # generate nonmissing genotype data 
        gen1, gen2, loc, n_per_loc = (i[keep_idx] for i in (geno1, geno2, loci_names, n_per_loci))
        @inbounds shared_loci[i] = length(keep_idx)
        @inbounds for (j, mthd) in enumerate(method)
            @inbounds relate_vecs[j][i] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples, inbreeding = inbreeding)
            boot_out = bootstrap_genos_nonmissing(gen1, gen2, loc, n_per_loc, allele_frequencies, method = mthd, iterations = iterations, inbreeding = inbreeding)
            @inbounds boot_means[j][i], boot_medians[j][i], boot_ses[j][i], boot_CI[j][i] = bootstrap_summary(boot_out, interval)
        end
        next!(p)
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

function relatedness_boot_all(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975), inbreeding::Bool) where F
    loci_names = Symbol.(loci(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_samples = length(samples(data))
    allele_frequencies = allele_freq(data)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    boot_means, boot_medians, boot_ses = map(i -> deepcopy(relate_vecs), 1:3)
    boot_CI = map(i -> Vector{Union{Missing,Tuple{Float64,Float64}}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    @inbounds Base.Threads.@threads for i in 1:length(sample_pairs)
        @inbounds sample_1 = sample_pairs[i][1]
        @inbounds sample_2 = sample_pairs[i][2]
        @inbounds geno1 = popdata_idx[(sample_1,)].genotype
        @inbounds geno2 = popdata_idx[(sample_2,)].genotype
        # get index for genotype appearing missing in at least one individual in the pair
        keep_idx = nonmissings(geno1, geno2)
        # generate nonmissing genotype data 
        gen1, gen2, loc, n_per_loc = (i[keep_idx] for i in (geno1, geno2, loci_names, n_per_loci))
        @inbounds shared_loci[i] = length(keep_idx)
        @inbounds for (j, mthd) in enumerate(method)
            @inbounds relate_vecs[j][i] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loci, n_samples = n_samples, inbreeding = inbreeding)
            boot_out = bootstrap_genos_all(geno1, geno2, loci_names, n_per_loci, allele_frequencies, method = mthd, iterations = iterations, inbreeding = inbreeding)
            @inbounds boot_means[j][i], boot_medians[j][i], boot_ses[j][i], boot_CI[j][i] = bootstrap_summary(boot_out, interval)
        end
        next!(p)
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


function relatedness_boot_all(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975), inbreeding::Bool) where F
    loci_names = Symbol.(loci(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_samples = length(samples(data))
    allele_frequencies = allele_freq(data)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
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

            # get index for genotype appearing missing in at least one individual in the pair
            keep_idx = nonmissings(geno1, geno2)
            # generate nonmissing genotype data 
            gen1, gen2, loc, n_per_loc = (i[keep_idx] for i in (geno1, geno2, loci_names, n_per_loci))

            # populate shared_loci array
            @inbounds shared_loci[idx] = length(loc)

            @inbounds for (i, mthd) in enumerate(method)
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loci, n_samples = n_samples, inbreeding = inbreeding)
                boot_out = bootstrap_genos_all(geno1, geno2, loci_names, n_per_loci, allele_frequencies, method = mthd, iterations = iterations, inbreeding = inbreeding)
                @inbounds boot_means[i][idx], boot_medians[i][idx], boot_ses[i][idx], boot_CI[i][idx] = bootstrap_summary(boot_out, interval)
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