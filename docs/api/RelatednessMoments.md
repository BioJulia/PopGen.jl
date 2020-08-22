---
id: relatednessmoments
title: RelatednessMoments.jl
sidebar_label: RelatednessMoments.jl
---

### `Blouin`
    Blouin(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Blouin (1996)
- Single Locus Equation: The number of alleles shared between individuals over ploidy.
    - If both allele positions are shared (e.g. AA x AA or AB x AB) then 1
    - If one allele position is shared (e.g. AB x AC) then 0.5
    - If neither allele position is shared (e.g. AB x CD) then 0
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Blouin, M. S., Parsons, M., Lacaille, V., & Lotz, S. (1996). Use of microsatellite loci to classify individuals by relatedness. Molecular ecology, 5(3), 393-401.


----

### `LiHorvitz`
    LiHorvitz(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Li and Horvitz (1953)
- Single Locus Equation: If all alleles are the same between individuals (eg. AA x AA) then 1.
    - If two alleles are shared between individuals (eg.  AA x AB or AB x AB) then 0.5.
    - If only one allele is shared between individuals (eg. AB x AC) then 0.25.
    - If no alleles are shared (eg. AB x CD) then 0.
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Li, C. C., & Horvitz, D. G. (1953). Some methods of estimating the inbreeding coefficient. American journal of human genetics, 5(2), 107.


----

### `Loiselle`
    Loiselle(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness using the estimator propsed by
Loiselle et al (1995) and modified to individual dyads by Heuertz et al. (2003).
- Multiple Locus Equation:
- Assumes no inbreeding

See equations 22 in: Wang(2017) for variant of estimator used

Loiselle, B. A., Sork, V. L., Nason, J., & Graham, C. (1995). Spatial genetic structure of a tropical understory shrub, _Psychotria officinalis_ (Rubiaceae). American journal of botany, 82(11), 1420-1425.

Heuertz, M., Vekemans, X., Hausman, J. F., Palada, M., & Hardy, O. J. (2003). Estimating seed vs. pollen dispersal from spatial genetic structure in the common ash. Molecular Ecology, 12(9), 2483-2495.

Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.


----

### `Lynch`
    Lynch(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Allele sharing index described by Lynch (1988)
- Single Locus Equation: If all alleles are the same between individuals (eg. AA x AA) then 1.
    - If both individuals are heterozygous with the same alleles or one is homozygous for the shared allele (eg. AB x AB or AA x AB) then 0.75.
    - If only one allele is shared between individuals (eg. AB x AC) then 0.5.
    - If no alleles are shared (eg. AB x CD) then 0.
- How to combine multiple loci: Single locus estimates are simply averaged together
- Assumes no inbreeding

Lynch, M. (1988). Estimation of relatedness by DNA fingerprinting. Molecular biology and evolution, 5(5), 584-599.


----

### `LynchLi`
    LynchLi(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Lynch (1988) & improved by Li et al. (1993).
- Single Locus Equation:
- How to combine multiple loci: Sum the difference between observed and expected similarity across all loci and then divide by the sum of 1 - the expected similarity
- Assumes no inbreeding

See equations 13 - 16 in Wang (2017) for variant of estimator used

Li, C. C., Weeks, D. E., & Chakravarti, A. (1993). Similarity of DNA fingerprints due to chance and relatedness. Human heredity, 43(1), 45-52.

Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.


----

### `LynchRitland`
    LynchRitland(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Lynch and Ritland (1999).
- Single Locus Equation:
- How to combine multiple loci: Weighted average of each term seperately weighted by the sample variance (assuming zero relatedness) and subsequently divided by the average sampling variance
- Assumes no inbreeding

See equation 10 in Wang (2017) for variant of estimator used

Lynch, M., & Ritland, K. (1999). Estimation of pairwise relatedness with molecular markers. Genetics, 152(4), 1753-1766.

Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.


----

### `Moran`
    Moran(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Reinterpretation of Moran's I (commonly used for spatial autocorrelation) to estimate genetic relatedness
by Hardy and Vekemans (1999)
- Multiple Locus Equation:
- Assumes no inbreeding

Hardy, O. J., & Vekemans, X. (1999). Isolation by distance in a continuous population: reconciliation between spatial autocorrelation analysis and population genetics models. Heredity, 83(2), 145-154.


----

### `QuellerGoodnight`
    QuellerGoodnight(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness developed by Queller & Goodnight (1989).
- Single Locus Equation:
- How to combine multiple loci:
    - Multiple loci are combined by independently summing the two numerator and two denominator terms before performing the final division and averaging the two components.
- Assumes no inbreeding

See equation 3 in Wang (2017) for variant of estimator used.

Queller, D. C., & Goodnight, K. F. (1989). Estimating relatedness using genetic markers. Evolution, 43(2), 258-275.

Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.


----

### `Ritland`
    Ritland(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness proposed by Li and Horvitz (1953) and implemented/made popular by Ritland (1996).
- Single Locus Equation:
- How to combine multiple loci: A weighted average of individual locus specific estimates weighted by sampling variance
- Assumes no inbreeding

See equation 7 in: Wang (2017) for variant of estimator used

Ritland, K. (1996). Estimators for pairwise relatedness and individual inbreeding coefficients. Genetics Research, 67(2), 175-185.

Wang, J. (2017). Estimating pairwise relatedness in a small sample of individuals. Heredity, 119(5), 302-313.


----

### `Wang`
    Wang(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the moments based estimator of pairwise relatedness by Wang (2002).
- Single Locus Equation:
- How to combine multiple loci: Each individual locus subcomponent (b-g) and each genotypic state (P1-P3) is averaged weighted by the average similarity of unrelated dyads at each locus. Then the values of V, Φ, Δ, and r are calculated
- Assumes no inbreeding
- Corrected for sampling bias in allele frequencies to get an unbiased estimator

Wang, J. (2002). An estimator for pairwise relatedness using molecular markers. Genetics, 160(3), 1203-1215.
