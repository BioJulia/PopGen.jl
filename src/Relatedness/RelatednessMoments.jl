export QuellerGoodnight, Ritland, Lynch, LynchRitland, LynchLi, LiHorvitz, Moran, Blouin, Loiselle, Wang, relatedness

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
function a_wang_base(m::Int, alleles::Dict)
    sum(values(alleles) .^ m)
end

function a_wang(N::Int, alleles::Dict)
    #unbiased estimator
    a = 0.0

    b = (N * a_wang_base(2, alleles) - 1) / (N - 1)

    c = (N^2 * a_wang_base(3, alleles) - 3 * (N - 1) * b - 1) / ((N - 1) * (N - 2))

    d = (N^3 * a_wang_base(4, alleles) - 6 * (N - 1) * (N - 2) * c - 7 * (N - 1) * b - 1) / (N^3 - 6 * N^2 + 11 * N - 6)

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

        #N = nonmissing(data.loci[data.loci.locus .== string(loc), :genotype])

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


"""
    relatedness(data::PopData; method::Function, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String)
    relatedness(data::PopData; method::Vector{Function}, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String)
    # to specify individuals for comparison
    relatedness(data::PopData, sample_names::Vector{String}; method::Function, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String)
    relatedness(data::PopData, sample_names::Vector{String}; method::Vector{Function}, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String)

Return a dataframe of pairwise relatedness estimates for all or select pairs of samples in a `PopData` object. **Note:** samples must be diploid.

### estimator methods
There are several estimators available and are listed below. `relatedness` takes the
function names as arguments (**case sensitive**), therefore do not use quotes or colons
in specifying the methods. Methods can be supplied as a vector.
- `Blouin`
- `LiHorvitz`
- `Loiselle`
- `Lynch`
- `LynchLi`
- `LynchRitland`
- `Moran`
- `QuellerGoodnight`
- `Ritland`
- `Wang`

For more information on a specific method, please see the respective docstring (e.g. `?Loiselle`).

### bootstrap methods
To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.05, 0.95)` (i.e. 90%),
however that can be changed by supplying a `Tuple` of `(low, high)` to the keyword `interval`.
There are two available resampling methods, `"all"` (default) and `"nonmissing"`.
- `"all"` : resamples all loci for a pair of individuals and then drops missing loci between them (recommended)
    - speed: slower
    - pro: better resampling variation
    - con: by chance some iterations may have a lot of missing loci that have to be dropped
- `"nonmissing"` : resamples only the shared non-missing loci between the pair
    - speed: faster
    - pro: every iteration guarantees the same number of loci compared between the pair
    - con: too-tight confidence intervals due to less possible variation

**Examples**
```
julia> cats = nancycats();

julia> relatedness(cats, method = Ritland, iterations = 100);
27966×8 DataFrame. Omitted printing of 2 columns
│ Row   │ sample_1 │ sample_2 │ n_loci │ Ritland    │ Ritland_mean │ Ritland_median │
│       │ String   │ String   │ Int64  │ Float64?   │ Float64?     │ Float64?       │
├───────┼──────────┼──────────┼────────┼────────────┼──────────────┼────────────────┤
│ 1     │ N215     │ N216     │ 8      │ 0.258824   │ 0.274537     │ 0.266457       │
│ 2     │ N215     │ N217     │ 8      │ 0.193238   │ 0.191591     │ 0.180821       │
│ 3     │ N215     │ N218     │ 8      │ 0.127497   │ 0.119988     │ 0.105559       │
│ 4     │ N215     │ N219     │ 8      │ 0.0453471  │ 0.0558557    │ 0.0573132      │
│ 5     │ N215     │ N220     │ 8      │ 0.108251   │ 0.12274      │ 0.110878       │
│ 6     │ N215     │ N221     │ 8      │ 0.205139   │ 0.206509     │ 0.204635       │
⋮
│ 27961 │ N297     │ N281     │ 7      │ -0.0487076 │ -0.052506    │ -0.0549532     │
│ 27962 │ N297     │ N289     │ 7      │ 0.154745   │ 0.170053     │ 0.169637       │
│ 27963 │ N297     │ N290     │ 7      │ 0.189647   │ 0.19302      │ 0.189994       │
│ 27964 │ N281     │ N289     │ 8      │ 0.0892068  │ 0.0914532    │ 0.0775897      │
│ 27965 │ N281     │ N290     │ 7      │ 0.104614   │ 0.107821     │ 0.115087       │
│ 27966 │ N289     │ N290     │ 7      │ 0.0511663  │ 0.0498539    │ 0.0547733      │

julia> relatedness(cats, ["N7", "N111", "N115"], method = [Ritland, Wang]);
3×5 DataFrame
│ Row │ sample_1 │ sample_2 │ n_loci │ Ritland    │ Wang      │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?  │
├─────┼──────────┼──────────┼────────┼────────────┼───────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.129432  │ -0.395806 │
│ 2   │ N7       │ N115     │ 9      │ -0.0183925 │ 0.0024775 │
│ 3   │ N111     │ N115     │ 9      │ 0.0240152  │ 0.183966  │

julia> relatedness(cats, ["N7", "N111", "N115"], method = [Loiselle, Moran], iterations = 100, interval = (0.025, 0.975));
3×13 DataFrame. Omitted printing of 7 columns
│ Row │ sample_1 │ sample_2 │ n_loci │ Loiselle   │ Loiselle_mean │ Loiselle_median │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?      │ Float64?        │
├─────┼──────────┼──────────┼────────┼────────────┼───────────────┼─────────────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.101618  │ 0.0141364     │ 0.00703954      │
│ 2   │ N7       │ N115     │ 9      │ -0.0428898 │ 0.0743497     │ 0.0798708       │
│ 3   │ N111     │ N115     │ 9      │ 0.13681    │ 0.266043      │ 0.257748        │
```
"""
function relatedness(data::PopData, sample_names::Vector{String}; method::F, iterations::Int64 = 0, interval::Tuple{Float64, Float64} = (0.025, 0.975), resample::String = "all") where F
    all(data.meta[data.meta.name .∈ Ref(sample_names), :ploidy] .== 2) == false && error("Relatedness analyses currently only support diploid samples")
    errs = ""
    all_samples = samples(data)
    if sample_names != all_samples
        for i in sample_names
            if i ∉ all_samples
                errs *= " $i,"
            end
        end
        errs != "" && error("Samples not found in the PopData: " * errs)
    end
    if eltype(method) != Function
        method = [method]
    end
    for i in Symbol.(method)
        if i ∉ [:QuellerGoodnight, :Ritland, :Lynch, :LynchLi, :LynchRitland, :Wang, :Loiselle, :Blouin, :Moran, :LiHorvitz, :dyadicLikelihood]
            errs *= "$i is not a valid method\n"
        end
    end
    errs != "" && error(errs * "Methods are case-sensitive. Please see the docstring (?relatedness) for additional help.")
    if iterations > 0
        if resample == "all"
            relatedness_boot_all(data, sample_names, method = method, iterations = iterations, interval = interval)
        elseif resample == "nonmissing"
            relatedness_boot_nonmissing(data, sample_names, method = method, iterations = iterations, interval = interval)
        else
            throw(ArgumentError("Please choose from resample methods \"all\" or \"nonmissing\""))
        end
    else
        relatedness_no_boot(data, sample_names, method = method)
    end
end


function relatedness(data::PopData; method::F, iterations::Int64 = 0, interval::Tuple{Float64, Float64} = (0.025, 0.975), resample::String = "all") where F
    sample_names = samples(data) |> collect
    relatedness(data, sample_names, method = method, iterations = iterations, interval = interval, resample = resample)
end