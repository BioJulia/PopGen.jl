"""
    bootstrap_summary(boot_out::Vector{Union{Missing, Float64}}, width::Tuple{Float64, Float64})

Return the mean, median, standard error, and quantiles (given by `witdth`) of relatedness resampling.
"""
@inline function bootstrap_summary(boot_out::Vector{Union{Missing, Float64}}, width::Tuple{Float64, Float64})
    all(ismissing.(boot_out)) == true && return missing, missing, missing, missing
    boot_skipmissing = collect(skipmissing(boot_out))
    n_nonmiss = length(boot_skipmissing)
    Mean = mean(boot_skipmissing)
    Median = median(boot_skipmissing)
    SE = sqrt(sum((boot_skipmissing - (boot_skipmissing / n_nonmiss)).^2) / (n_nonmiss - 1))
    quants = quantile(boot_skipmissing, width)

    return Mean, Median, SE, quants
end

"""
    bootstrap_genos_all(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::NamedTuple; method::Function, iterations::Int)

Perform `iterations` number of bootstrap resampling iterations of all genotypes between pair (`ind1` `ind2`). Returns a vector of length `interatotions`
of the relatedness estimate given by method `method`. This is an internal function with `locus_names`, `n_per_loc`, and `alleles` supplied by `relatedness_boot_all`.
"""
@inline function bootstrap_genos_all(ind1::T, ind2::T, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::U; method::F, iterations::Int) where T <: GenoArray where U <: NamedTuple where F
    relate_vec_boot = Vector{Union{Missing,Float64}}(undef, iterations)
    n_loc = length(locus_names)
    for iter in 1:iterations
        # bootstrap the indices
        boot_idx = rand(1:n_loc, n_loc)
        # sample the source vectors with the resampled/bootstrapped indices
        ind1_boot, ind2_boot, loc_boot, n_per_loci = map(i -> getindex(i, boot_idx), [ind1, ind2, locus_names, n_per_loc])
        # get index for genotype appearing missing in at least one individual in the pair
        keep_idx = nonmissings(ind1_boot, ind2_boot)
        relate_vec_boot[iter] = method(ind1_boot[keep_idx], ind2_boot[keep_idx], loc_boot[keep_idx], alleles, loc_n = n_per_loci[keep_idx], n_samples = n_loc)
    end
    return relate_vec_boot
end


"""
    bootstrap_genos_nonmissing(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::NamedTuple; method::Function, iterations::Int)

Perform `iterations` number of bootstrap resampling iterations of only shared (nonmissing) genotypes between pair (`ind1` `ind2`). Returns a vector of length `interatotions`
of the relatedness estimate given by method `method`. This is an internal function with `locus_names`, `n_per_loc`, and `alleles` supplied by `relatedness_boot_nonmissing`.
"""
@inline function bootstrap_genos_nonmissing(ind1::T, ind2::T, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::U; method::F, iterations::Int) where T <: GenoArray where U <: NamedTuple where F
    relate_vec_boot = Vector{Union{Missing,Float64}}(undef, iterations)
    n_loc = length(locus_names)
    for iter in 1:iterations
        # bootstrap the indices
        boot_idx = rand(1:n_loc, n_loc)
        # sample the source vectors with the resampled/bootstrapped indices
        ind1_boot, ind2_boot, loc_boot, n_per_loci = map(i -> getindex(i, boot_idx), [ind1, ind2, locus_names, n_per_loc])
        # faster/cheaper n counting
        relate_vec_boot[iter] = method(ind1_boot, ind2_boot, loc_boot, alleles, loc_n = n_per_loci, n_samples = n_loc)
    end
    return relate_vec_boot
end

"""
    relatedness_boot_all(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64}) where F

Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. Bootstrapping resamples using
the `all` method, where resampling occurs over all loci. This is an internal function with all arguments provided by `relatedness`.
"""
function relatedness_boot_all(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975)) where F
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
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loci, n_samples = n_samples)
                boot_out = bootstrap_genos_all(geno1, geno2, loci_names, n_per_loci, allele_frequencies, method = mthd, iterations = iterations)
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


"""
    relatedness_boot_nonmissing(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64}) where F

Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. Bootstrapping resamples using
the `nonmissing` method, where resampling occurs over only shared non-missing loci. This is an internal function with all arguments provided by `relatedness`.
"""
function relatedness_boot_nonmissing(data::PopData, sample_names::Vector{String}; method::F, iterations::Int = 100, interval::Tuple{Float64, Float64} = (0.025, 0.975)) where F
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
        @inbounds geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            @inbounds geno2 = popdata_idx[(ind2,)].genotype

            # get index for genotype appearing missing in at least one individual in the pair
            keep_idx = nonmissings(geno1, geno2)
            # generate nonmissing genotype data 
            gen1, gen2, loc, n_per_loc = (i[keep_idx] for i in (geno1, geno2, loci_names, n_per_loci))

            # populate shared_loci array
            @inbounds shared_loci[idx] = length(loc)
            
            @inbounds for (i, mthd) in enumerate(method)
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples)
                boot_out = bootstrap_genos_nonmissing(gen1, gen2, loc, n_per_loc, allele_frequencies, method = mthd, iterations = iterations)
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

"""
    relatedness_no_boot(data::PopData, sample_names::Vector{String}; method::F) where F

Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. 
This is an internal function with arguments provided by `relatedness`.
"""
function relatedness_no_boot(data::PopData, sample_names::Vector{String}; method::F) where F
    loci_names = Symbol.(loci(data))
    n_samples = length(samples(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        @inbounds geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            @inbounds geno2 = popdata_idx[(ind2,)].genotype

            # get index for genotype appearing missing in at least one individual in the pair
            keep_idx = nonmissings(geno1, geno2)
            
            # populate shared_loci array
            @inbounds shared_loci[idx] = length(keep_idx)
            @inbounds [relate_vecs[i][idx] = @inbounds mth(geno1[keep_idx], geno2[keep_idx], loci_names[keep_idx], allele_frequencies, loc_n = n_per_loci[keep_idx], n_samples = n_samples) for (i,mth) in enumerate(method)]

            update!(p, idx)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci)
    [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    return out_df
end