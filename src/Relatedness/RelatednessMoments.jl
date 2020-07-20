#### Simple code versions ####
"""
    QuellerGoodnight(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness developed by Queller & Goodnight (1989).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 3 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function QuellerGoodnight(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple

    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    numerator2 = 0.0
    denominator1 = 0.0
    denominator2 = 0.0

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
        a,b = gen1
        c,d = gen2
        ident = ((a == c) + (a == d) + (b == c) + (b == d))
        fq_a = alleles[loc][a]
        fq_b = alleles[loc][b]
        fq_c = alleles[loc][c]
        fq_d = alleles[loc][d]
        numerator1 += ident - 2.0 * (fq_a + fq_b)
        numerator2 += ident - 2.0 * (fq_c + fq_d)

        denominator1 += (2.0 * (1.0 + (a==b) - fq_a - fq_b))
        denominator2 += (2.0 * (1.0 + (c==d) - fq_c - fq_d))
    end
    return (numerator1/denominator1 + numerator2/denominator2)/2.0
end

"""
    Ritland(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness proposed by Li and Horvitz (1953) and implemented/made popular by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 7 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function Ritland(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple

    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    denominator1 = 0.0
    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
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
    LynchRitland(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness by Ritland (1996).
- Bases allele frequencies on entire population
- Inbreeding can only be assumed not to exist.
See equation 10 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
Ritland original citation: https://www.cambridge.org/core/journals/genetics-research/article/estimators-for-pairwise-relatedness-and-individual-inbreeding-coefficients/9AE218BF6BF09CCCE18121AA63561CF7
"""
function LynchRitland(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    denominator1 = 0.0

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
        a,b = gen1
        c,d = gen2
        fq_a = alleles[loc][a]
        fq_b = alleles[loc][b]
        fq_c = alleles[loc][c]
        fq_d = alleles[loc][d]

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
    LynchLi(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness by Lynch (1988) & improved by Li et al. (1993).
See equations 13 - 16 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function LynchLi(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    denominator1 = 0.0

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
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
    Loiselle(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness using the estimator propsed by
Loiselle et al (1995) and modified to individual dyads by Heuertz et al. (2003).
See equations 22 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function Lioselle(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    numerator1 = 0.0
    denominator1 = 0.0
    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
        a,b = gen1
        c,d = gen2

        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            n += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            d += fq * (1.0 - fq)
        end
    end
    return (2.0 * numerator1 / denominator1) + 2.0 / (2 * length(samples(data)) - 1)
end

### Wang 2002 helper functions ###
function a_wang(m::Int, alleles::Dict)::Float64
    #TODO Change to unbiased formulation (eq 25)
    sum(values(alleles) .^ m)
end

function b_wang(alleles::Dict)
    2.0 * a_wang(2, alleles) ^ 2.0 - a_wang(4, alleles)
end

function c_wang(alleles::Dict)
    a_wang(2, alleles) - 2.0 * a_wang(2, alleles) ^ 2
end

function d_wang(alleles::Dict)
    4.0 * (a_wang(3, alleles) - a_wang(4, alleles))
end

function e_wang(alleles::Dict)
    2.0 * (a_wang(2, alleles) - 3.0 * a_wang(3, alleles) + 2.0 * a_wang(4, alleles))
end

function f_wang(alleles::Dict)
    4.0 * (a_wang(2, alleles) - a_wang(2, alleles)^2 - 2.0 * a_wang(3, alleles) + 2.0 * a_wang(4, alleles))
end

function g_wang(alleles::Dict)
    1.0 - 7.0 * a_wang(2, alleles) + 4.0 * a_wang(2, alleles)^2 + 10.0 * a_wang(3, alleles) - 8.0 * a_wang(4, alleles)
end

function u_wang(alleles::Dict)
    2.0 * a_wang(2, alleles) - a_wang(3, alleles)
end

"""
    Wang(data::PopData, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness by Wang (2002).
See https://www.genetics.org/content/genetics/160/3/1203.full.pdf
"""
function Wang(data::PopData, ind1::String, ind2::String; alleles::T) where T <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS

    r = Vector{Float64}(undef, length(loci(data)))
    u = Vector{Float64}(undef, length(loci(data)))

    geno1 = get_genotypes(data, ind1)
    geno2 = get_genotypes(data, ind2)

    loc_id = 0
    for (loc,gen1,gen2) in zip(skipmissings(Symbol.(loci(data)), geno1, geno2)...)
        loc_id += 1
        i,j = gen1
        k,l = gen2

        b = b_wang(alleles[loc])
        c = c_wang(alleles[loc])
        d = d_wang(alleles[loc])
        e = e_wang(alleles[loc])
        f = f_wang(alleles[loc])
        g = g_wang(alleles[loc])
        u[loc_id] = u_wang(alleles[loc])

        # Which category of dyad
        # Both alleles shared between individuals either the same or different
        P1 = 1.0 * ((i == j == k == l) | (i == k & j == l) | (i == l & k == j))
        # One allele shared between individuals and one is homozygous for that allele
        P2 = 1.0 * ((i == j == k != l) | (i == j == l != k) | (k == l == i != j) | (k == l == j != i))
        # One allele shared with the other two being unique
        P3 = 1.0 * ((i == k + i != j + i != l + j != l) | (i == l + i != k + i != j + j != k) | (j == k + j != i + j != l + l != i) | (j == l + j != i + j != k + k != l))
        P4 = 1.0 * ((P1 + P2 + P3) == 0)

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

        #Eq 1.0
        r[loc_id] = (Φ/2.0 + Δ)
    end
    return 1 / (sum(1/u) * u) * r
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
    relate_vec = zeros(length(sample_pairs))
    idx = 0
    @inbounds for (sample_n, ind1) in enumerate(sample_names[1:end-1])
        #@inbounds for ind2 in sample_names[sample_n+1:end]
        #@inbounds Base.Threads.@threads for ind2 in sample_names[sample_n+1:end]
        @inbounds @sync Base.Threads.@spawn for ind2 in sample_names[sample_n+1:end]
            idx += 1
            #TODO Add column for number of shared and number of missing loci for each pair
            relate_vec[idx] += method(data, ind1, ind2, alleles = allele_frequencies)
        end
    end
    method_colname = Symbol("$method")
    return DataFrame(:sample_1 => getindex.(sample_pairs, 1), :sample_2 => getindex.(sample_pairs, 2), method_colname => relate_vec)
end
