## THIS IS A SCRATCH-WORK FILE  -- PAVEL ##

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



function relatedness_no_boot_fix(data::PopData, sample_names::Vector{String}; method::F, inbreeding::Bool) where F
    loci_names = Symbol.(loci(data))
    n_samples = length(samples(data))
    sample_pairs = pairwise_pairs(sample_names)
    n_per_loci = DataFrames.combine(groupby(data.loci, :locus), :genotype => nonmissing => :n)[:, :n]
    allele_frequencies = allele_freq(data)
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    p = Progress(length(sample_pairs), dt = 1, color = :blue)
    popdata_idx = groupby(data.loci, :name)
    idx = Base.Threads.Atomic{Int}()
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        @inbounds geno1 = popdata_idx[(ind1,)].genotype
        @inbounds @sync Base.Threads.@threads for ind2 in sample_names[sample_n+1:end]
            Base.Threads.atomic_add!(idx, 1)
            @inbounds geno2 = popdata_idx[(ind2,)].genotype

            # get index for genotype appearing missing in at least one individual in the pair
            keep_idx = nonmissings(geno1, geno2)
            
            # populate shared_loci array
            @inbounds shared_loci[idx.value] = length(keep_idx)
            @inbounds [relate_vecs[i][idx.value] = @inbounds mth(geno1[keep_idx], geno2[keep_idx], loci_names[keep_idx], allele_frequencies, loc_n = n_per_loci[keep_idx], n_samples = n_samples, inbreeding = inbreeding) for (i,mth) in enumerate(method)]

            update!(p, idx.value)
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci)
    [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    return out_df
end


function relatedness2(data::PopData, sample_names::Vector{String}; method::F, iterations::Int64 = 0, interval::Tuple{Float64, Float64} = (0.025, 0.975), resample::String = "all", inbreeding::Bool = false) where F
    all(data.meta[data.meta.name .∈ Ref(sample_names), :ploidy] .== 2) == false && error("Relatedness analyses currently only support diploid samples")
    errs = ""
    all_samples = samples(data)
    if sample_names != all_samples
        [errs *= "$i," for i in sample_names if i ∉ all_samples]
        errs != "" && error("Samples not found in the PopData: " * errs)
    end
    if eltype(method) != Function
        method = [method]
    end
    relate_mthds = [:QuellerGoodnight, :Ritland, :Lynch, :LynchLi, :LynchRitland, :Wang, :Loiselle, :Blouin, :Moran, :LiHorvitz, :dyadicLikelihood]
    [errs *= "$i is not a valid method\n" for i in Symbol.(method) if i ∉ relate_mthds]
    errs != "" && throw(ArgumentError(errs * "Methods are case-sensitive. Please see the docstring (?relatedness) for additional help."))
    if iterations > 0
        if resample == "all"
            relatedness_boot_all(data, sample_names, method = method, iterations = iterations, interval = interval, inbreeding = inbreeding)
        elseif resample == "nonmissing"
            relatedness_boot_nonmissing(data, sample_names, method = method, iterations = iterations, interval = interval, inbreeding = inbreeding)
        else
            throw(ArgumentError("Please choose from resample methods \"all\" or \"nonmissing\""))
        end
    else
        relatedness_no_boot_fix2(data, sample_names, method = method, inbreeding = inbreeding)
    end
end


function relatedness2(data::PopData; method::F, iterations::Int64 = 0, interval::Tuple{Float64, Float64} = (0.025, 0.975), resample::String = "all", inbreeding::Bool = false) where F
    sample_names = samples(data) |> collect
    relatedness2(data, sample_names, method = method, iterations = iterations, interval = interval, resample = resample, inbreeding = inbreeding)
end


## original function -- threading does not work as expected ###
function relatedness_no_boot(data::PopData, sample_names::Vector{String}; method::F, inbreeding::Bool) where F
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
            @inbounds [relate_vecs[i][idx] = @inbounds mth(geno1[keep_idx], geno2[keep_idx], loci_names[keep_idx], allele_frequencies, loc_n = n_per_loci[keep_idx], n_samples = n_samples, inbreeding = inbreeding) for (i,mth) in enumerate(method)]

            update!(p, idx)
        end
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
                @inbounds relate_vecs[i][idx] = mthd(gen1, gen2, loc, allele_frequencies, loc_n = n_per_loc, n_samples = n_samples, inbreeding = inbreeding)
                boot_out = bootstrap_genos_nonmissing(gen1, gen2, loc, n_per_loc, allele_frequencies, method = mthd, iterations = iterations, inbreeding = inbreeding)
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
















# method using kwargs and Symbols as methods. Incomplete.
function relatedness(data::PopData, sample_names::Vector{String}; kwargs...)
    #method::F, iterations::Int64 = 0, interval::Tuple{Float64, Float64} = (0.025, 0.975), resample::String = "all") where F
    kw_dict = Dict(kwargs...)
    all(data.meta[data.meta.name .∈ Ref(sample_names), :ploidy] .== 2) == false && error("Relatedness analyses currently only support diploid samples")
    all_samples = samples(data)
    errs = ""
    if sample_names != all_samples
        [errs *= "$i," for i in sample_names if i ∉ all_samples]
        errs != "" && error("Samples not found in the PopData: " * errs)
    end
    method = kw_dict[:method]
    if !(typeof(method) <: AbstractVector)
        method = [method]
    end
    !(haskey(kw_dict, :method)) && throw(ArgumentError("Please supply the keyword argument: method. See the docstring (?relatedness) for additional help"))
    relate_mthds = [:QuellerGoodnight, :Ritland, :Lynch, :LynchLi, :LynchRitland, :Wang, :Loiselle, :Blouin, :Moran, :LiHorvitz, :dyadicLikelihood]
    [errs *= "$i is not a valid method\n" for i in method if i ∉ relate_mthds]
    errs != "" && throw(ArgumentError(errs * "Methods are case-sensitive. Please see the docstring (?relatedness) for additional help."))
    if haskey(kw_dict, :iterations)
        if kw_dict[:iterations] > 0
            if haskey(kw_dict, :resample)
                "all" != kw_dict[:resample] != "nonmissing" && throw(ArgumentError("Please choose from resample methods \"all\" or \"nonmissing\""))
                kw_dict[:resample] == "all" && relatedness_boot_all(data, sample_names, method = method, iterations = iterations, interval = interval)
                kw_dict[:resample] == "nonmissing" && relatedness_boot_nonmissing(data, sample_names, method = method, iterations = iterations, interval = interval)
            else
                relatedness_boot_all(data, sample_names, method = method, iterations = iterations, interval = interval)
            end
    else
        relatedness_no_boot(data, sample_names, method = method)
    end
end



#= functions as Symbols for arguments

julia> function do_thing(x, operation::Symbol)
           getfield(@__MODULE__, operation)(x)
       end

do_thing (generic function with 1 method)

julia> do_thing(1, :sin)
0.8414709848078965

julia> do_thing(1, :sdfgsdfsdf)
ERROR: UndefVarError: sdfgsdfsdf not defined
Stacktrace:
 [1] do_thing(::Int64, ::Symbol) at ./REPL[26]:2
 [2] top-level scope at REPL[28]:1
=#





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