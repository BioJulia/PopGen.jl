---
id: relatedness
title: Relatedness (Kinship)
sidebar_label: Relatedness (Kinship)
---

**Background**

Sometimes you want or need to know the relatedness of individuals in a sample. Relatedness is exactly what its name implies: how related individuals are given some provided genetic information (e.g. full siblings, half siblings, etc.). Relatedness can be used in quantitative genetics to estimate heritability, additive genetic variances, and covariances. It can also be used in population genetics to study isolation-by-distance or population structure.

The goal of calculating relatedness from molecular markers is to accurately estimate the proportion of the genome which is identical by descent between two individuals. With a pedigree this is "relatively" straightforward. However, for large, natural, populations pedigrees tend not to exist and so many brilliant minds have developed various equations and algorithms to estimate the relatedness from a set of molecular markers. Given two diploid individuals, there are 9 "identity by descent" models available between them ([Jacquard 1975](https://www.springer.com/gp/book/9783642884177), paywall), as shown below (from [Milligan 2003](https://www.genetics.org/content/163/3/1153.full)):

![Jacquard IBD](/PopGen.jl/img/jacquard_identitiies.jpg)

Broadly speaking there are two different ways of estimating genetic relatedness using molecular markers, methods of moments, and likelihood estimators. Generally moments estimators will be faster but aren't constrained to being between the theoretical minimum and maximum values of 0 and 1. The likelihood estimators use likelihood functions derived from the set of Jacquard Identity States (above) to determine the most likely inheritance pattern. One difference between the two classes is [generally] moments estimators require an assumption of no inbreeding, while that assumption isn't necessarily required for likelihood estimators (though it does simplify the math). It is increasingly common to use multiple estimators on pairs, simulated from your molecular marker, with a known relationships to determine the most appropriate estimator to use with your given data.

PopGen.jl currently implements one of each class of estimator. The moments based estimator developed by [Queller & Goodnight 1989](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x) using the variant defined by [Lynch & Ritland 1999](https://www.genetics.org/content/152/4/1753.full) and the likelihood estimator created by [Milligan 2003](https://www.genetics.org/content/163/3/1153.full). You can imagine there's a lot that  happens under the hood to perform this for all loci across all individuals-- all of which dutifully written by Jason Selwyn :star2:.

:::note A note about removing kin
There are reasons for removing kin in population genetics datasets. For one, there are no siblings/kin or mixed-generations in a Hardy-Weinberg Equilibrium population, and the inclusion of siblings/kin in analyses that rely on HWE assumptions [technically] violate those assumptions. However, there are also arguments to keep kin/siblings in your data, those data are important for effective population size, linkage disequilibrium, etc. 



see  [Waples and Anderson (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14022)
::: 

## Calculate Relatedness

Calculate the relatedness of all pairs of individuals in the dataframe using either Milligan's Dyadic Maximum Likelihood estimator or Queller & Goodnight's estimator.

```julia
pairwise_relatedness(data::PopObj; method::String, inbreeding::Bool = true, verbose::Bool = true)
```

### arguments

- `data` : the input `PopData`

### keyword arguments

- `method` : Method of relatedness estimation (see below)
- `inbreeding` : Include the possibility of inbreeding (true) or not (false) - Only used with `method = "dyadml"`
- `verbose` : If `false` only progress bar will be shown. If `true` extra output will be shown depending on the method chosen

### methods

- `"dyadml"` : Milligan (2003) Dyadic Likelihood Relatedness
- `"qg"` : Queller & Goodnight (1989) Relatedness

*****

### example
:::: tabs card stretch
::: tab Queller-Goodnight relatedness
```julia
cats = nancycats() ;
pairwise_relatedness(cats, method = "qg", verbose = false)
```
:::
::: tab output
```

```

See [Development](/other/api/hidden_api.md) for all of the APIs associated with `relatedness()`

## Speed

Depending on the number samples and loci in your data, the maximum-likelihood approach to relatedness can be quite time consuming. We include a progress bar thanks to [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) to provide some indication of how long it will take. As a point of reference, it takes approximately 13 hours to perform this relatedness calculation on the `gulfsharks` data (212 samples x 2213 loci).

Currently, relatedness calculations run single-threaded, and we hope to parallelize it with the stable release of Julia 1.3  to make it even faster. Many hands make light work!

-------------------

## Acknowledgements

Both [Convex.jl](https://github.com/JuliaOpt/Convex.jl) and [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl) are pivotal for these calculations, and we thank the authors for their time developing and maintaining them, along with the members of the Julia Slack channel for pointing us towards those packages.
