---
id: relatedness
title: Relatedness (Kinship)
sidebar_label: Relatedness (Kinship)
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />

## Background

Sometimes you want or need to know the relatedness of individuals in a sample. Relatedness is exactly what its name implies: how related individuals are given some provided genetic information (e.g. full siblings, half siblings, etc.). Relatedness can be used in quantitative genetics to estimate heritability, additive genetic variances, and covariances. It can also be used in population genetics to study isolation-by-distance or population structure.

The goal of calculating relatedness from molecular markers is to accurately estimate the proportion of the genome which is identical by descent between two individuals. With a pedigree this is "relatively" straightforward. However, for large, natural, populations pedigrees tend not to exist and some brilliant minds have developed various equations to estimate the relatedness from a set of molecular markers. Given two diploid individuals, there are 9 "identity by descent" models available between them ([Jacquard 1975](https://www.springer.com/gp/book/9783642884177), paywall), as shown below (from [Milligan 2003](https://www.genetics.org/content/163/3/1153.full)):

![Jacquard IBD](/PopGen.jl/img/jacquard_identitiies.jpg)

Broadly speaking there are two different ways of estimating genetic relatedness using molecular markers, methods of moments, and likelihood estimators. Generally, moments estimators will be faster but aren't constrained to being between the theoretical minimum and maximum values of 0 and 1. The likelihood estimators use likelihood functions derived from the set of Jacquard Identity States (above) to determine the most likely inheritance pattern. One difference between the two classes is [generally] moments estimators require an assumption of no inbreeding, while that assumption isn't necessarily required for likelihood estimators (though it does simplify the math). It is increasingly common to use multiple estimators on pairs, simulated from your molecular markers, with a known relationships to determine the most appropriate estimator to use with your given data.

PopGen.jl implements a wide variety of moments-based estimators: Blouin, Li & Horvitz, Loiselle, Lynch, Lynch/Li, Lynch & Ritland, Moran, Queller & Goodnight, Ritland, and Wang. Along with these, we provide an option to estimate mean, median, standard error, and confidence intervals using bootstrapping.

:::note A note about removing kin
There are reasons for removing kin in population genetics datasets. For one, there are no siblings/kin or mixed-generations in a Hardy-Weinberg Equilibrium population, and the inclusion of siblings/kin in analyses that rely on HWE assumptions [technically] violate those assumptions. However, there are also arguments to keep kin/siblings in your data, those data are important for effective population size, linkage disequilibrium, etc. 

see  [Waples and Anderson (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14022)
::: 

## Calculate Relatedness
<Tabs
  block={true}
  defaultValue="all"
  values={[
    { label: 'All vs. All', value: 'a', },
    { label: 'Specific Samples', value: 's', },
  ]
}>
<TabItem value="a">
```julia
relatedness(data::PopData; method::Function, iterations::Int64, interval::Tuple(Float64, Float64))
relatedness(data::PopData; method::Vector{Function}, iterations::Int64, interval::Tuple(Float64, Float64))
```
Return a dataframe of pairwise relatedness estimates for all individuals in a `PopData` object.
To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.275, 0.975)` (i.e. 95%), however that can be changed by supplying a `Tuple{Float64, Float64}` of `(low, high)` to the keyword `interval`. **Note:** samples must be diploid.

</TabItem>
<TabItem value="s">
```julia
relatedness(data::PopData, sample_names::Vector{String}; method::Function, iterations::Int64, interval::Tuple(Float64, Float64))
relatedness(data::PopData, sample_names::Vector{String}; method::Vector{Function}, iterations::Int64, interval::Tuple(Float64, Float64))
```
Return a dataframe of pairwise relatedness estimates for all pairs of the supplied sample names in a `PopData` object.
To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.275, 0.975)` (i.e. 95%),
however that can be changed by supplying a `Tuple{Float64, Float64}` of `(low, high)` to the keyword `interval`.
**Note:** samples must be diploid.

</TabItem>
</Tabs>

**Methods**

There are several estimators available and are listed below. `relatedness` takes the
function names as arguments (**case sensitive**), therefore do not use quotes or colons
in specifying the methods. Methods can be supplied as a vector. 

:::tip autocompletion
Since the methods correspond to function names, they will tab-autocomplete when 
inputting them. For more information on a specific method, please see the respective docstring (e.g. `?
Loiselle`).
:::

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

**Examples**

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'Single Method', value: 's', },
    { label: 'Multiple Methods', value: 'm', },
    { label: 'With Bootstrapping', value: 'b', },
  ]
}>
<TabItem value="s">

```julia
julia> cats = nancycats();

julia> cat_kin = relatendess(cats, samples(cats)[1:10], method = Ritland)
```

</TabItem>
<TabItem value="m">

```julia
julia> cats = nancycats();

julia> cat_kin = relatendess(cats, samples(cats)[1:10], method = [Moran, QuellerGoodnight])
```

</TabItem>
<TabItem value="b">

```julia
julia> cats = nancycats();

julia> cat_kin = relatendess(cats, samples(cats)[1:10], method = Ritland, iterations = 100)
```

</TabItem>
</Tabs>

Read **Relatedness Estimators** below for more in-depth information regarding each estimator.

## Relatedness Estimators
### Blouin
The moments based estimator developed by [Blouin (year)](). Call `method = Blouin` to use this method. 

### Li & Horvitz
The moments based estimator developed by [Li & Horvitz (year)](). Call `method = LiHorvitz` to use this method. 

### Loiselle
The moments based estimator developed by [Loiselle (year)](). Call `method = Loiselle` to use this method. 

### Lynch
The moments based estimator developed by [Lynch (year)](). Call `method = Lynch` to use this method. 

### Lynch / Li
The moments based estimator developed by [Lynch/Li (year)](). Call `method = LynchLi` to use this method. 

### Lynch & Ritland
The moments based estimator developed by [Lynch & Ritland (year)](). Call `method = LynchRitland` to use this method. 

### Moran
The moments based estimator developed by [Moran (year)](). Call `method = Moran` to use this method. 

### Queller & Goodnight
The moments based estimator developed by [Queller & Goodnight (1989)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x). Call `method = QuellerGoodnight` to use this method. 


### Ritland
The moments based estimator developed by [Ritland (year)](). Call `method = Ritland` to use this method. 

### Wang
The moments based estimator developed by [Wang (year)](). Call `method = Wang` to use this method. 


---------------------
## Acknowledgements

