"""
    Blouin(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Blouin (1996)

- Single Locus Equation: The number of alleles shared between individuals over ploidy.
    - If both allele positions are shared (e.g. AA x AA or AB x AB) then 1
    - If one allele position is shared (e.g. AB x AC) then 0.5
    - If neither allele position is shared (e.g. AB x CD) then 0
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Blouin, M. S., Parsons, M., Lacaille, V., & Lotz, S. (1996). Use of microsatellite loci to classify individuals by relatedness. Molecular ecology, 5(3), 393-401.
"""
function Blouin(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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
    LiHorvitz(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Li and Horvitz (1953)

-Single Locus Equation: If all alleles are the same between individuals (eg. AA x AA) then 1.
    - If two alleles are shared between individuals (eg.  AA x AB or AB x AB) then 0.5.
    - If only one allele is shared between individuals (eg. AB x AC) then 0.25.
    - If no alleles are shared (eg. AB x CD) then 0.
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Li, C. C., & Horvitz, D. G. (1953). Some methods of estimating the inbreeding coefficient. American journal of human genetics, 5(2), 107.
"""
function LiHorvitz(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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
    Loiselle(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness using the estimator propsed by
Loiselle et al (1995) and modified to individual dyads by Heuertz et al. (2003).

- Multiple Locus Equation:
- Assumes no inbreeding

See equations 22 in: Wang(2017) for variant of estimator used

Loiselle, B. A., Sork, V. L., Nason, J., & Graham, C. (1995). Spatial genetic structure of a tropical understory shrub, <i>Psychotria officinalis</i> (Rubiaceae). American journal of botany, 82(11), 1420-1425.
Heuertz, M., Vekemans, X., Hausman, J. F., Palada, M., & Hardy, O. J. (2003). Estimating seed vs. pollen dispersal from spatial genetic structure in the common ash. Molecular Ecology, 12(9), 2483-2495.
Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.
"""
function Loiselle(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    isempty(locus_names) && return missing
    d_kw = Dict(kwargs...)
    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            numerator1 += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            denominator1 += fq * (1.0 - fq)
        end
    end
    return numerator1 / denominator1 + 2.0 / (2 * d_kw[:n_samples] - 1)
end


"""
    Lynch(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Lynch (1988)

- Single Locus Equation: If all alleles are the same between individuals (eg. AA x AA) then 1.
    - If both individuals are heterozygous with the same alleles or one is homozygous for the shared allele (eg. AB x AB or AA x AB) then 0.75.
    - If only one allele is shared between individuals (eg. AB x AC) then 0.5.
    - If no alleles are shared (eg. AB x CD) then 0.
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Lynch, M. (1988). Estimation of relatedness by DNA fingerprinting. Molecular biology and evolution, 5(5), 584-599.
"""
function Lynch(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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
    LynchLi(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Lynch (1988) & improved by Li et al. (1993).

- Single Locus Equation:
- How to combine multiple loci: Sum the difference between observed and expected similarity across all loci and then divide by the sum of 1 - the expected similarity
- Assumes no inbreeding

See equations 13 - 16 in Wang (2017) for variant of estimator used

Li, C. C., Weeks, D. E., & Chakravarti, A. (1993). Similarity of DNA fingerprints due to chance and relatedness. Human heredity, 43(1), 45-52.
Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.
"""
function LynchLi(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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
    LynchRitland(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Lynch and Ritland (1999).

- Single Locus Equation:
- How to combine multiple loci: Weighted average of each term seperately weighted by the sample variance (assuming zero relatedness) and subsequently divided by the average sampling variance
- Assumes no inbreeding

See equation 10 in Wang (2017) for variant of estimator used

Lynch, M., & Ritland, K. (1999). Estimation of pairwise relatedness with molecular markers. Genetics, 152(4), 1753-1766.
Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.
"""
function LynchRitland(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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
    Moran(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Reinterpretation of Moran's I (commonly used for spatial autocorrelation) to estimate genetic relatedness
by Hardy and Vekemans (1999)

- Multiple Locus Equation:
- Assumes no inbreeding

Hardy, O. J., & Vekemans, X. (1999). Isolation by distance in a continuous population: reconciliation between spatial autocorrelation analysis and population genetics models. Heredity, 83(2), 145-154.
"""
function Moran(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    isempty(locus_names) && return missing

    numerator1 = 0.0
    denominator1 = 0.0

    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            numerator1 += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            #denominator1 += ((sum(gen1 .== allele) / 2.0) - fq)^2
            denominator1 += (((sum(gen1 .== allele) / 2.0) - fq)^2 + ((sum(gen2 .== allele) / 2.0) - fq)^2) / 2.0
        end
        #denominator1 += (1 / (length(alleles[loc]) - 1))
    end
    return (numerator1 / denominator1)
end

function Moran_experimental(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    isempty(locus_names) && return missing

    numerator1 = Vector{Float64}(undef, length(locus_names))
    denominator1 = similar(numerator1)

    numerator1 = numerator1 .* 0.0
    denominator1 = denominator1 .* 0.0

    idx = 0
    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        idx += 1
        for allele in keys(alleles[loc])
            fq = alleles[loc][allele]
            numerator1[idx] += ((sum(gen1 .== allele) / 2.0) - fq) * ((sum(gen2 .== allele) / 2.0) - fq)
            denominator1[idx] += (((sum(gen1 .== allele) / 2.0) - fq)^2) #+ ((sum(gen2 .== allele) / 2.0) - fq)^2)  / 2
        end
        denominator1[idx] += (1 / (length(alleles[loc]) - 1))
    end

    return mean(numerator1 ./ denominator1)
end

"""
    QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness developed by Queller & Goodnight (1989).

- Single Locus Equation:
- How to combine multiple loci:
    - Multiple loci are combined by independently summing the two numerator and two denominator terms before performing the final division and averaging the two components.
- Assumes no inbreeding
See equation 3 in Wang(2017) for variant of estimator used.

Queller, D. C., & Goodnight, K. F. (1989). Estimating relatedness using genetic markers. Evolution, 43(2), 258-275.
Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.
"""
function QuellerGoodnight(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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

- Single Locus Equation:
- How to combine multiple loci: A weighted average of individual locus specific estimates weighted by sampling variance
- Assumes no inbreeding

See equation 7 in: Wang (2017) for variant of estimator used

Ritland, K. (1996). Estimators for pairwise relatedness and individual inbreeding coefficients. Genetics Research, 67(2), 175-185.
Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.
"""
function Ritland(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
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

### Wang 2002 helper functions ###
function _a_wang_base(m::Int, alleles::Dict)
    sum(values(alleles) .^ m)
end

function _a_wang(N::Int, alleles::Dict)
    #unbiased estimator
    a = 0.0

    b = (N * _a_wang_base(2, alleles) - 1) / (N - 1)

    c = (N^2 * _a_wang_base(3, alleles) - 3 * (N - 1) * b - 1) / ((N - 1) * (N - 2))

    d = (N^3 * _a_wang_base(4, alleles) - 6 * (N - 1) * (N - 2) * c - 7 * (N - 1) * b - 1) / (N^3 - 6 * N^2 + 11 * N - 6)

    return [a, b, c, d]
end


"""
Wang(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Wang (2002).

-Single Locus Equation:
-How to combine multiple loci: Each individual locus subcomponent (b-g) and each genotypic state (P1-P3) is averaged weighted by the average similarity of unrelated dyads at each locus. Then the values of V, Φ, Δ, and r are calculated

-Assumes no inbreeding
-Corrected for sampling bias in allele frequencies to get an unbiased estimator

Wang, J. (2002). An estimator for pairwise relatedness using molecular markers. Genetics, 160(3), 1203-1215.
"""
function Wang(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    #TODO NEED TO CHECK TO CONFIRM EQUATIONS
    isempty(locus_names) && return missing
    kw_dict = Dict(kwargs...)
    P1 = Vector{Float64}(undef, length(locus_names))
    P2, P3, P4, u, b, c, d, e, f, g = map(i -> similar(P1), 1:10)
    loc_id = 0

    for (loc,gen1,gen2, N) in zip(locus_names, ind1, ind2, kw_dict[:loc_n])
        loc_id += 1
        i,j = gen1
        k,l = gen2

        #N = nonmissing(data.genodata[data.genodata.locus .== string(loc), :genotype])

        a = _a_wang(2 * N, alleles[loc])
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
