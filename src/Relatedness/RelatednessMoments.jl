#### Simple code versions ####
"""
    QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Calculates the moments based estimator of pairwise relatedness developed by Queller & Goodnight (1989).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 3 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function QuellerGoodnight(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    numerator1 = 0.0
    numerator2 = 0.0
    denominator1 = 0.0
    denominator2 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        a,b = gen1
        c,d = gen2
        ident = ((a == c) + (a == d) + (b == c) + (b == d))
        fq_a, fq_b, fq_c, fq_d = map(i -> alleles[loc][i], (a,b,c,d))

        numerator1 += ident - 2.0 * (fq_a + fq_b)
        numerator2 += ident - 2.0 * (fq_c + fq_d)

        denominator1 += (2.0 * (1.0 + (a==b) - fq_a - fq_b))
        denominator2 += (2.0 * (1.0 + (c==d) - fq_c - fq_d))
    end
    return (numerator1/denominator1 + numerator2/denominator2)/2.0
end

"""
    Ritland(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Calculates the moments based estimator of pairwise relatedness proposed by Li and Horvitz (1953) and implemented/made popular by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 7 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function Ritland(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        a,b = gen1
        c,d = gen2

        A = ((alleles[loc] |> length) - 1)

        R = 0.0
        for i in unique((a,b,c,d))
            # Individual locus relatedness value (eq 7 in paper)
            R += ((((a == i) + (b == i)) * ((c == i) + (d == i))) / (4.0 * alleles[loc][i]))
        end
        R = (2.0 / A) * (R - 1.0)
        # numerator for weighted combination of loci
        numerator1 += (R * A)
        # denominator for weighted combination of loci
        denominator1 += A
    end
    return numerator1 / denominator1
end

"""
    LynchRitland(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Calculates the moments based estimator of pairwise relatedness by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 10 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function LynchRitland(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        a,b = gen1
        c,d = gen2
        fq_a, fq_b, fq_c, fq_d = map(i -> alleles[loc][i], (a,b,c,d))

        n1 = fq_a * ((b == c) + (b == d)) + fq_b * ((a == c) + (a == d)) - 4.0 * fq_a * fq_b
        n2 = fq_c * ((d == a) + (d == b)) + fq_d * ((c == a) + (c == b)) - 4.0 * fq_c * fq_d

        d1 = 2.0 * (1.0 + (a == b)) * (fq_a + fq_b) - 8.0 * fq_a * fq_b
        d2 = 2.0 * (1.0 + (c == d)) * (fq_c + fq_d) - 8.0 * fq_c * fq_d


        WL1 = ((1 + (a == b)) * (fq_a + fq_b) - 4 * fq_a * fq_b) / (2 * fq_a * fq_b)
        WL2 = ((1 + (c == d)) * (fq_c + fq_d) - 4 * fq_c * fq_d) / (2 * fq_c * fq_d)

        numerator1 += ((n1 / d1) * WL1 + (n2 / d2) * WL2)
        denominator1 += (WL1 + WL2) / 2.0
    end
    return numerator1 / denominator1
end

"""
    LynchLi(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Calculates the moments based estimator of pairwise relatedness by Lynch (1988) & improved by Li et al. (1993).
See equations 13 - 16 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function LynchLi(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        a,b = gen1
        c,d = gen2

        Sxy = (0.5) * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) + ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
        #TODO Change to unbiased formulation (eq 25)
        S0 = 2.0 * sum(values(alleles[loc]) .^ 2) - sum(values(alleles[loc]) .^ 3)

        numerator1 += Sxy - S0
        denominator1 += 1.0 - S0
    end
    return numerator1 / denominator1
end

"""
    Loiselle(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Calculates the moments based estimator of pairwise relatedness using the estimator propsed by
Loiselle et al (1995) and modified to individual dyads by Heuertz et al. (2003).
See equations 22 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function Loiselle(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            numerator1 += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            denominator1 += fq * (1.0 - fq)
        end
    end
    return numerator1 / denominator1 + 2.0 / (2 * length(samples(data)) - 1)
end

"""
    LiHorvitz(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Allele sharing index described by Li and Horvitz (1953)
"""
function LiHorvitz(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    Bxy = Vector{Float64}(undef, length(locus_names))

    loc_id = 0
    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        loc_id += 1
        i,j = gen1
        k,l = gen2

        Bxy[loc_id] = sum([i, j] .∈ [k,l]') / 4
    end
    return mean(Bxy)
end

"""
    Lynch(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Allele sharing index described by Lynch (1988)
"""
function Lynch(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    Sxy = Vector{Float64}(undef, length(locus_names))
    loc_id = 0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        loc_id += 1
        i,j = gen1
        k,l = gen2

        Sxy[loc_id] = ((i ∈ gen2) + (j ∈ gen2) + (k ∈ gen1) + (l ∈ gen1)) / 4
    end
    return mean(Sxy)
end

"""
    Blouin(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Allele sharing index described by Blouin (1996)
"""
function Blouin(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing

    Mxy = Vector{Float64}(undef, length(locus_names))
    loc_id = 0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        loc_id += 1
        i,j = gen1
        k,l = gen2

        Mxy[loc_id] = (((i ∈ gen2) & (k ∈ gen1)) + ((j ∈ gen2) & (l ∈ gen1))) / 2
    end
    return mean(Mxy)
end

"""
Moran(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple) 
Reinterpretation of Moran's I (commonly used for spatial autocorrelation) to estimate genetic relatedness
by Hardy and Vekemans (1999)
"""
function Moran(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            numerator1 += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            denominator1 += ((sum(gen1 .== allele) / 2.0) - fq)^2
        end
        denominator1 += (1 / (length(alleles[loc]) - 1))
    end
    return (numerator1 / denominator1)
end

### Wang 2002 helper functions ###
function a_wang_base(m::Int, alleles::Dict)
    sum(values(alleles) .^ m)
end

function a_wang(N::Int, alleles::Dict)
    #unbiased estimator
    a = Vector{Float64}(undef, 4)
    a[1] = 0.0

    a[2] = (N * a_wang_base(2, alleles) - 1) / (N - 1)

    a[3] = (N^2 * a_wang_base(3, alleles) - 3 * (N - 1) * a[2] - 1) / ((N - 1) * (N - 2))

    a[4] = (N^3 * a_wang_base(4, alleles) - 6 * (N - 1) * (N - 2) * a[3] - 7 * (N - 1) * a[2] - 1) / (N^3 - 6 * N^2 + 11 * N - 6)

    return a
end

"""
Wang(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Wang (2002).
See https://www.genetics.org/content/genetics/160/3/1203.full.pdf
"""
function Wang(ind1::T, ind2::T, locus_names::Vector{Symbol}; alleles::U) where T <: GenoArray where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    isempty(locus_names) && return missing
    P1 = Vector{Float64}(undef, length(loci(data)))
    P2, P3, P4, u, b, c, d, e, f, g = map(i -> similar(P1), 1:10)
    loc_id = 0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        loc_id += 1
        i,j = gen1
        k,l = gen2

        N = nonmissing(data.loci[data.loci.locus .== string(loc), :genotype])

        a = a_wang(2 * N, alleles[loc])
        a2_sq = a[2] ^ 2

        u[loc_id] = 2 * a[2] - a[3]

        # Which category of dyad
        Sxy = ((i ∈ gen2) + (j ∈ gen2) + (k ∈ gen1) + (l ∈ gen1)) / 4

        # Both alleles shared between individuals either the same or different
        P1[loc_id] = 1.0 * (Sxy == 1)
        # One allele shared between individuals and one is homozygous for that allele
        P2[loc_id] = 1.0 * (Sxy == (3/4))
        # One allele shared with the other two being unique
        P3[loc_id] = 1.0 * (Sxy == (1/2))
        P4[loc_id] = 1.0 * ((P1 + P2 + P3) == 0)

        b[loc_id] = (2.0 * a2_sq - a[4])
        c[loc_id] = (a[2] - 2.0 * a2_sq + a[4])
        d[loc_id] = (4.0 * (a[3] - a[4]))
        e[loc_id] = (2.0 * (a[2] - 3.0 * a[3] + 2.0 * a[4]))
        f[loc_id] = (4.0 * (a[2] - a2_sq - 2.0 * a[3] + 2.0 * a[4]))
        g[loc_id] = (1.0 - 7.0 * a[2] + 4.0 * a2_sq + 10.0 * a[3] - 8.0 * a[4])

    end
    #return (1 / (sum(1/u) * u)) * r
    w = (1 / (sum(1/u) * u))


    P1 = w * P1
    P2 = w * P2
    P3 = w * P3

    b = w * b
    c = w * c
    d = w * d
    e = w * e
    f = w * f
    g = w * g

    #Eq 11
    V = (1.0 - b)^2 * (e^2 * f + d * g^2) -
        (1.0 - b) * (e * f - d * g)^2 +
        2.0 * c * d * f * (1.0 - b) * (g + e) +
        c^2 * d * f * (d + f)

    #Eq 9
    Φ = (d * f * ((e + g) * (1.0 - b) + c * (d + f)) * (P1 - 1.0) +
        d * (1.0 - b) * (g * (1.0 - b - d) + f * (c + e)) * P3 +
        f * (1.0 - b) * (e * (1.0 - b - f) + d * (c + g)) * P2) / V

    #Eq 10
    Δ = (c * d * f * (e + g) * (P1 + 1.0 - 2 * b) +
        ((1.0 - b) * (f * e^2 + d * g^2) - (e * f - d * g)^2) * (P1 - b) +
        c * (d * g - e * f) * (d * P3 - f * P2) - c^2 * d * f * (P3 + P2 - d - f) -
        c * (1.0 - b) * (d * g * P3 + e * f * P2)) / V

    r = (Φ/2.0 + Δ)
    return r
end


#TODO this is 100% incomplete
function pairwise_relatedness(data::PopData; method::Function, inbreeding::Bool = true, verbose::Bool = true)
    # check that dataset is entirely diploid
    all(data.meta.ploidy .== 2) == false && error("Relatedness analyses currently only support diploid samples")

    allele_frequencies = NamedTuple{Tuple(Symbol.(loci(data)))}(
                            Tuple(allele_freq.(locus.(Ref(data), loci(data))))
                        )
    sample_names = samples(data)
    sample_pairs = [tuple(sample_names[i], sample_names[j]) for i in 1:length(sample_names)-1 for j in i+1:length(sample_names)]
    relate_vec = Vector{Union{Missing,Float64}}(undef, length(sample_pairs))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = get_genotypes(data, ind1)
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            #TODO Progress Bar
            geno2 = get_genotypes(data, ind2)

            # filter out loci missing in at least one individual in the pair
            loc,gen1,gen2 = collect.(skipmissings(Symbol.(loci(data)), geno1, geno2))

            # populate shared_loci array
            shared_loci[idx] = length(loc)
            relate_vec[idx] = method(gen1, gen2, loc, alleles = allele_frequencies)

            #TODO Bootstrap loop
            #TODO Bootstrap post-process
        end
    end
    method_colname = Symbol("$method")
    return DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :shared_loci => shared_loci, method_colname => relate_vec)
end



function pairwise_relatedness(data::PopData; method::Union{Function, Vector{Function}}, inbreeding::Bool = true, verbose::Bool = true)
    # check that dataset is entirely diploid
    all(data.meta.ploidy .== 2) == false && error("Relatedness analyses currently only support diploid samples")

    allele_frequencies = NamedTuple{Tuple(Symbol.(loci(data)))}(
                            Tuple(allele_freq.(locus.(Ref(data), loci(data))))
                        )
    sample_names = samples(data)
    sample_pairs = [tuple(sample_names[i], sample_names[j]) for i in 1:length(sample_names)-1 for j in i+1:length(sample_names)]
    
    if eltype(method) != Function
        method = [method]
    end
    relate_vecs = map(i -> Vector{Union{Missing,Float64}}(undef, length(sample_pairs)), 1:length(method))
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = get_genotypes(data, ind1)
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            #TODO Progress Bar
            geno2 = get_genotypes(data, ind2)

            # filter out loci missing in at least one individual in the pair
            loc,gen1,gen2 = collect.(skipmissings(Symbol.(loci(data)), geno1, geno2))

            # populate shared_loci array
            shared_loci[idx] = length(loc)
            [relate_vecs[i][idx] = mth(gen1, gen2, loc, alleles = allele_frequencies) for (i,mth) in enumerate(method)]

            #TODO Bootstrap loop
            #TODO Bootstrap post-process
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci)
    [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    return out_df
end


#dictionary implementation
function pairwise_relatedness3(data::PopData; method::Vector{Function}, inbreeding::Bool = true, verbose::Bool = true)
    # check that dataset is entirely diploid
    all(data.meta.ploidy .== 2) == false && error("Relatedness analyses currently only support diploid samples")

    allele_frequencies = NamedTuple{Tuple(Symbol.(loci(data)))}(
                            Tuple(allele_freq.(locus.(Ref(data), loci(data))))
                        )
    sample_names = samples(data)
    sample_pairs = [tuple(sample_names[i], sample_names[j]) for i in 1:length(sample_names)-1 for j in i+1:length(sample_names)]
    relate_vecs = Dict([i => Vector{Union{Missing,Float64}}(undef, length(sample_pairs)) for i in Symbol.(method)]...)
    shared_loci = Vector{Int}(undef, length(sample_pairs))
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        geno1 = get_genotypes(data, ind1)
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            #TODO Progress Bar
            geno2 = get_genotypes(data, ind2)

            # filter out loci missing in at least one individual in the pair
            loc,gen1,gen2 = collect.(skipmissings(Symbol.(loci(data)), geno1, geno2))

            # populate shared_loci array
            shared_loci[idx] = length(loc)
            [relate_vecs[Symbol(mthd)][idx] = mthd(gen1, gen2, loc, alleles = allele_frequencies) for mthd in method]

            #TODO Bootstrap loop
            #TODO Bootstrap post-process
        end
    end
    method_colnames = [Symbol("$i") for i in method]
    out_df = DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), :n_loci => shared_loci, relate_vecs...)
   # [out_df[:, mth] = relate_vecs[i] for (i, mth) in enumerate(method_colnames)]
    #return out_df
end