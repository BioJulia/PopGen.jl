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

![Jacquard IBD](/img/jacquard_identitiies.jpg)

Broadly speaking there are two different ways of estimating genetic relatedness using molecular markers: methods of moments, and likelihood estimators. Generally, moments estimators will be faster but aren't constrained to being between the theoretical minimum and maximum values of 0 and 1. The likelihood estimators use likelihood functions derived from the set of Jacquard Identity States (above) to determine the most likely inheritance pattern. One difference between the two classes is (generally) moments estimators require an assumption of no inbreeding, while that assumption isn't necessarily required for likelihood estimators (though it does simplify the math). It is increasingly common to use multiple estimators on pairs, simulated from your molecular markers, with known relationships to determine the most appropriate estimator to use with your data.

PopGen.jl implements a wide variety of moments-based estimators: Blouin, Li & Horvitz, Loiselle, Lynch, Lynch/Li, Lynch & Ritland, Moran, Queller & Goodnight, Ritland, and Wang. Along with these, we provide an option to estimate mean, median, standard error, and confidence intervals using bootstrapping. Have a look at [the guide](/blog/relatedness) we provide detailing how to perform a basic relatedness analysis.

:::note A note about removing kin
There are reasons for removing kin in population genetics datasets. For one, there are no siblings/kin or mixed-generations in a Hardy-Weinberg Equilibrium population, and the inclusion of siblings/kin in analyses that rely on HWE assumptions [technically] violate those assumptions. However, there are also arguments to keep kin/siblings in your data, those data are important for effective population size, linkage disequilibrium, etc. 

see  [Waples and Anderson (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14022)
::: 

## Calculate Relatedness

<Tabs
  block={true}
  defaultValue="a"
  values={[
    { label: 'All vs. All', value: 'a', },
    { label: 'Specific Samples', value: 's', },
  ]
}>
<TabItem value="a">

```julia
relatedness(data::PopData; method::F, iterations::Int64, interval::Tuple(Float64, Float64))
```
Return a dataframe of pairwise relatedness estimates for all individuals in a `PopData` object using 
method `F` (see below). To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.275, 0.975)` (i.e. 95%), however that can be changed by supplying a `Tuple{Float64, Float64}` of `(low, high)` to the keyword `interval`. **Note:** samples must be diploid.

### Arguments
- `data` : A PopData object
- `sample_names` : A list of samples names to calculate relatedness for (optional)

### Keyword Arguments
- `method` : A method function or vector of method functions (see below)
- `iterations` : The number of iterations to perform bootstrapping (default: `0`, will not perform bootstrapping)
- `interval` : A Tuple of (low, high) indicating the confidence intervals you would like for bootstrapping (default: `(0.275, 0.975)`, i.e. 95%)
- `inbreeding` : true/false of whether to consider inbreeding in the calculations (default: `false`). Only used in `dyadML`

**Examples**
```
julia> cats = @nancycats;

julia> relatedness(cats, method = Ritland)
27966×4 DataFrame
    Row │ sample_1  sample_2  n_loci  Ritland
        │ String    String    Int64   Float64?
 ───────┼─────────────────────────────────────────
      1 │ N215      N216           8   0.258824
      2 │ N215      N217           8   0.193238
      3 │ N215      N218           8   0.127497
      4 │ N215      N219           8   0.0453471
   ⋮   │    ⋮         ⋮        ⋮          ⋮
  27963 │ N297      N290           7   0.189647
  27964 │ N281      N289           8   0.0892068
  27965 │ N281      N290           7   0.104614
  27966 │ N289      N290           7   0.0511663
                                27958 rows omitted
```

</TabItem>
<TabItem value="s">

```julia
relatedness(data::PopData, sample_names::Vector{String}; method::F, iterations::Int64, interval::Tuple(Float64, Float64))
```
Return a dataframe of pairwise relatedness estimates for all pairs of the supplied sample names in a `PopData` object using 
method `F` (see below). To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.275, 0.975)` (i.e. 95%),
however that can be changed by supplying a `Tuple{Float64, Float64}` of `(low, high)` to the keyword `interval`.
**Note:** samples must be diploid.

### Arguments
- `data` : A PopData object
- `sample_names` : A list of samples names to calculate relatedness for (optional)

### Keyword Arguments
- `method` : A method function or vector of method functions (see below)
- `iterations` : The number of iterations to perform bootstrapping (default: `0`, will not perform bootstrapping)
- `interval` : A Tuple of (low, high) indicating the confidence intervals you would like for bootstrapping (default: `(0.275, 0.975)`, i.e. 95%)
- `inbreeding` : true/false of whether to consider inbreeding in the calculations (default: `false`). Only used in `dyadML`

**Examples**
```
julia> cats = @nancycats;

julia> relatedness(cats, ["N7", "N111", "N115"], method = [Ritland, Wang])
3×5 DataFrame
│ Row │ sample_1 │ sample_2 │ n_loci │ Ritland    │ Wang      │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?  │
├─────┼──────────┼──────────┼────────┼────────────┼───────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.129432  │ -0.395806 │
│ 2   │ N7       │ N115     │ 9      │ -0.0183925 │ 0.0024775 │
│ 3   │ N111     │ N115     │ 9      │ 0.0240152  │ 0.183966  │


julia> relatedness(cats, ["N7", "N111", "N115"], method = [Loiselle, Moran], iterations = 100, interval = (0.025, 0.975))
3×13 DataFrame. Omitted printing of 7 columns
│ Row │ sample_1 │ sample_2 │ n_loci │ Loiselle   │ Loiselle_mean │ Loiselle_median │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?      │ Float64?        │
├─────┼──────────┼──────────┼────────┼────────────┼───────────────┼─────────────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.101618  │ 0.0141364     │ 0.00703954      │
│ 2   │ N7       │ N115     │ 9      │ -0.0428898 │ 0.0743497     │ 0.0798708       │
│ 3   │ N111     │ N115     │ 9      │ 0.13681    │ 0.266043      │ 0.257748        │
```

</TabItem>
</Tabs>

### Methods
There are several estimators available and are listed below. `relatedness` takes the
function names as arguments (**case sensitive**), therefore do not use quotes or colons
in specifying the methods. Methods can be supplied as a vector. 

- [Blouin](analyses/relatedness.md#blouin)
- [LiHorvitz](analyses/relatedness.md#li--horvitz)
- [Loiselle](analyses/relatedness.md#loiselle)
- [Lynch](analyses/relatedness.md#lynch)
- [LynchLi](analyses/relatedness.md#lynch--li)
- [LynchRitland](analyses/relatedness.md#lynch--ritland)
- [Moran](analyses/relatedness.md#lynch--moran)
- [QuellerGoodnight](analyses/relatedness.md#queller--goodnight)
- [Ritland](analyses/relatedness.md#ritland)
##### In Progress (incomplete*)
- [dyadML](analyses/relatedness.md#dyadic-maximum-likelihood)
- [Wang](analyses/relatedness.md#wang)

:::note *more relatedness
Contact us or submit a pull request if you're interested in contributing to the relatedness methods. We're currently in the process of adding dryadML and Wang (2002) estimators and speeding up the existing methods. More help is always welcomed! :smile:
:::

#### Examples

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
julia> cats = @nancycats;

julia> cat_kin = relatendess(cats, samples(cats)[1:10], method = Ritland)
```

</TabItem>
<TabItem value="m">

```julia
julia> cats = @nancycats;

julia> cat_kin = relatendess(cats, samples(cats)[1:10], method = [Moran, QuellerGoodnight])
```

</TabItem>
<TabItem value="b">

```julia
julia> cats = @nancycats;

julia> cat_kin = relatendess(cats, method = [Ritland, Wang], iterations = 100)
```

</TabItem>
</Tabs>

:::tip autocompletion
Since the methods correspond to function names, they will tab-autocomplete when 
inputting them. For more information on a specific method, please see the respective docstring (e.g. `?Loiselle`).
:::

## Relatedness Estimators
### Blouin
The moments based estimator developed by [Blouin (1996)](https://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.1996.00094.x). Call `method = Blouin` to use this method. 

### Li & Horvitz
The moments based estimator developed by [Li & Horvitz (1953)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1716461/). Call `method = LiHorvitz` to use this method. 

### Loiselle
The moments based estimator developed by [Loiselle (1995)](https://bsapubs.onlinelibrary.wiley.com/doi/abs/10.1002/j.1537-2197.1995.tb12679.x). Call `method = Loiselle` to use this method. 

### Lynch
The moments based estimator developed by [Lynch (1988)](https://pubmed.ncbi.nlm.nih.gov/3193879/). Call `method = Lynch` to use this method. 

### Lynch / Li
The moments based estimator developed by [Lynch/Li (1993)](https://pubmed.ncbi.nlm.nih.gov/8514326/). Call `method = LynchLi` to use this method. 

### Lynch & Ritland
The moments based estimator developed by [Lynch & Ritland (1999)](https://www.genetics.org/content/152/4/1753.short). Call `method = LynchRitland` to use this method. 

### Moran
The moments based estimator developed by [Moran (1950)](https://www.jstor.org/stable/2332142?origin=crossref&seq=1#metadata_info_tab_contents). Call `method = Moran` to use this method. 

### Queller & Goodnight
The moments based estimator developed by [Queller & Goodnight (1989)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x). Call `method = QuellerGoodnight` to use this method. 


### Ritland
The moments based estimator developed by [Ritland (1996)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1996.tb02347.x). Call `method = Ritland` to use this method. 

### Incomplete
#### Dyadic Maximum Likelihood
The moments based estimator developed by [Milligan (2003)](https://pubmed.ncbi.nlm.nih.gov/12663552/).

#### Wang
The moments based estimator developed by [Wang (2002)](https://www.genetics.org/content/160/3/1203.short).


## Posthoc analyses
There are several different kinds of things you can do with kinship information (e.g. network analysis), and one that's provided is lovingly
called `relatedness_posthoc()`, which performs a permutation analysis to
test if within-population relatedness is significantly greater than between-population relatedness for each population. We recommend that you
correct for multiple testing using `MultipleTesting.jl`.

```julia
relatedness_posthoc(data::PopData, results::Union{DataFrame, NamedTuple}; iterations::Int)
```
### Arguments
- `data` : A PopData object
- `results` : the DataFrame or NamedTuple results from `relatedness()`

### Keyword Arguments
- `iterations` : number of iterations for the permutation tests (default: `20000`)

:::tip not a great name
We admit "relatedness_posthoc" is not a great name for this function. Please
contact us with better ideas! :grin:
:::

**Example**
```
julia> cats = @nancycats ;

julia> rel_out = relatedness(cats, method = [Ritland, Moran], iterations = 100);

julia> relatedness_posthoc(cats, rel_out)
17x3 DataFrame
 Row │ population  Ritland_P  Moran_P
     │ String      Float64    Float64
─────┼────────────────────────────────
   1 │ 1              5.0e-5   5.0e-5
   2 │ 2              5.0e-5   5.0e-5
   3 │ 3              5.0e-5   5.0e-5
   4 │ 4              5.0e-5   5.0e-5
   5 │ 5              5.0e-5   5.0e-5
   6 │ 6              5.0e-5   5.0e-5
   7 │ 7              5.0e-5   5.0e-5
   8 │ 8              5.0e-5   5.0e-5
   9 │ 9              5.0e-5   5.0e-5
  10 │ 10             5.0e-5   5.0e-5
  11 │ 11             5.0e-5   5.0e-5
  12 │ 12             5.0e-5   5.0e-5
  13 │ 13             5.0e-5   5.0e-5
  14 │ 14             5.0e-5   5.0e-5
  15 │ 15             5.0e-5   5.0e-5
  16 │ 16             5.0e-5   5.0e-5
  17 │ 17             5.0e-5   5.0e-5
```

---------------------
## Acknowledgements
The relatedness methods were dutifully written and verified against R analogues by Jason Selwyn. These anaylses can take a while, especially if performing bootstrapping; we provide a progress bar via `ProgressMeter.jl` so you can move on and focus on other things in the meantime. 