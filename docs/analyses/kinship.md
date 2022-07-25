---
id: kinship
title: Kinship
sidebar_label: Kinship
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
There are two main methods of calculating pairwise relatedness, an all x all comparison of an entire PopData object or an all x all comparison
for a subset of individuals in a PopData object. Regardless of which you prefer, they can be perfomed without bootstrapping, returning a
`NamedMatrix`, or with bootstrapping, returning a `DataFrame`. **Neither method estimates self-relatedness, so the diagonals of the NamedMatrix should be ignored.** The resulting `NamedMatrix` can be converted to a `DataFrame` using `kinshiptotable()`.

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
kinship(data::PopData; kwargs...)
```
Calculate pairwise relatedness estimates for all individuals in a `PopData` object using 
the specified `method` (see below). Returns a `NamedMatrix` if not performing bootstrapping, otherwise returns a `DataFrame` (since bootstrapping provides more output information). To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `[0.0275, 0.975]` (i.e. 95%), however that can be changed by supplying a `Vector{Float64}` of `[low, high]` to the keyword `interval`. **Note:** samples must be diploid.

#### Arguments
- `data` : A PopData object

#### Keyword Arguments
- `method::Function` : A method function (see below)
- `iterations::Int64` : The number of iterations to perform bootstrapping (default: `0`, will not perform bootstrapping)
- `interval::Vector{Float64}` : A Vector of [low, high] indicating the confidence intervals you would like for bootstrapping (default: `[0.275, 0.975]`, i.e. 95%)

<Tabs
  block={true}
  defaultValue="wo"
  values={[
    { label: 'Without Bootstrapping', value: 'wo', },
    { label: 'With Bootstrapping', value: 'b', },
  ]
}>
<TabItem value="wo">

```julia
julia> cats = @nancycats ; 

julia> kin = kinship(cats, method = Moran)
237×237 Named Matrix{Float64}
A ╲ B │         N215          N216  …          N289          N290
──────┼──────────────────────────────────────────────────────────
N215  │ 8.13724e-316       1.62338  …       1.04589       1.15351
N216  │      1.62338       0.29485         0.957724        1.1637
N217  │     0.673577      0.587163         0.547427      0.709887
N218  │     0.896935       0.72942         0.919448      0.791255
⋮                  ⋮             ⋮  ⋱             ⋮             ⋮
N297  │     0.757915      0.858834          1.15432        1.2677
N281  │     0.686057      0.604236         0.942749       1.08762
N289  │      1.04589      0.957724              0.0         1.104
N290  │      1.15351        1.1637  …         1.104           0.0
```

</TabItem>
<TabItem value="b">

```julia
julia> cats = @nancycats ; 

julia> kinship_new(cats, method = Moran, iterations = 100)
27966×7 DataFrame
   Row │ sample1  sample2  Moran     bootmean  std       CI_lower      CI_upper 
       │ String   String   Float64   Float64   Float64   Float64       Float64  
───────┼────────────────────────────────────────────────────────────────────────
     1 │ N215     N216     1.62338   0.376626  0.27286    0.00274863   0.916719
     2 │ N215     N217     0.673577  0.202888  0.20094    0.00105976   0.59871
     3 │ N215     N218     0.896935  0.206272  0.232048   7.58373e-5   0.786113
     4 │ N215     N219     0.988931  0.236503  0.221345  -0.00053018   0.718204
   ⋮   │    ⋮        ⋮        ⋮         ⋮         ⋮           ⋮           ⋮
 27964 │ N281     N289     0.942749  0.220475  0.200358   0.001656     0.799307
 27965 │ N281     N290     1.08762   0.285053  0.289967   0.000299019  1.09343
 27966 │ N289     N290     1.104     0.277445  0.235519   0.00186445   0.858206
                                                              27959 rows omitted
```

</TabItem>
</Tabs>

</TabItem>
<TabItem value="s">

```julia
kinship(data::PopData, sample_names::Vector{AbstractString}; kwargs...)
```
Calculate pairwise relatedness estimates for all pairs of the supplied sample names in a `PopData` object using 
the specified `method` (see below). Returns a `NamedMatrix` if not performing bootstrapping, otherwise returns a `DataFrame` (since bootstrapping provides more output information). To calculate means, median, standard error, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `[0.0275, 0.975]` (i.e. 95%),
however that can be changed by supplying a `Vector{Float64}` of `[low, high]` to the keyword `interval`.
**Note:** samples must be diploid.

#### Arguments
- `data` : A PopData object
- `sample_names` : A list of samples names to calculate relatedness for

#### Keyword Arguments
- `method::Function` : A method function or vector of method functions (see below)
- `iterations::Int64` : The number of iterations to perform bootstrapping (default: `0`, will not perform bootstrapping)
- `interval::Vector{Float64}` : A Vector of [low, high] indicating the confidence intervals you would like for bootstrapping (default: `[0.0275, 0.975]`, i.e. 95%)


<Tabs
  block={true}
  defaultValue="wo"
  values={[
    { label: 'Without Bootstrapping', value: 'wo', },
    { label: 'With Bootstrapping', value: 'b', },
  ]
}>
<TabItem value="wo">

```julia
julia> cats = @nancycats ; 

julia> kin = kinship(cats, samplenames(cats)[1:10], method = Moran) ;
10×10 Named Matrix{Float64}
A ╲ B │         N215          N216  …          N223          N224
──────┼──────────────────────────────────────────────────────────
N215  │     5.0e-323       2.40275  …       1.28205      0.793651
N216  │      2.40275      5.0e-324          1.16525       1.03383
N217  │     0.860832      0.755287         0.938453      0.543519
N218  │      1.17754      0.967118         0.824974       1.00656
N219  │      1.36268       1.80995         0.760537       0.96706
N220  │      1.47059       1.84275         0.735041       1.39631
N221  │      2.03837       1.96335          1.11537       1.14702
N222  │     0.657895      0.804829          0.78841      0.905063
N223  │      1.28205       1.16525         5.0e-324      0.654129
N224  │     0.793651       1.03383  …      0.654129  6.93947e-310
```

</TabItem>
<TabItem value="b">

```julia
julia> cats = @nancycats ; 

julia> kin = kinship(cats, samplenames(cats)[1:10], method = Moran, iteratons = 100)
45×7 DataFrame
 Row │ sample1  sample2  Moran     bootmean   std       CI_lower      CI_upper 
     │ String   String   Float64   Float64    Float64   Float64       Float64  
─────┼─────────────────────────────────────────────────────────────────────────
   1 │ N215     N216     2.40275   0.166804   0.192498  -0.00103402   0.558154
   2 │ N215     N217     0.860832  0.132997   0.157386   0.00079658   0.500806
   3 │ N215     N218     1.17754   0.138991   0.193997   0.00178298   0.620476
   4 │ N215     N219     1.36268   0.118936   0.169152   0.000739799  0.527514
  ⋮  │    ⋮        ⋮        ⋮          ⋮         ⋮           ⋮           ⋮
  43 │ N222     N223     0.78841   0.0947265  0.1087    -0.000671851  0.316755
  44 │ N222     N224     0.905063  0.121569   0.125794   0.000570054  0.387836
  45 │ N223     N224     0.654129  0.0729529  0.121482   0.00105724   0.403654
                                                                38 rows omitted
```

</TabItem>
</Tabs>

</TabItem>
</Tabs>

### Methods
There are several estimators available and are listed below. `kinship` takes the
function names as arguments (**case sensitive**), therefore do not use quotes or colons
in specifying the methods.

:::tip autocompletion
Since the methods correspond to function names, they will tab-autocomplete when 
inputting them. For more information on a specific method, please see the respective docstring (e.g. `?Loiselle`).
:::

| Method | Type | Method Call |
|:----|:-----|:-----|
| [Blouin 1996](https://onlinelibrary.wiley.com/doi/10.1046/j.1365-294X.1996.00094.x) | moments-based | `Blouin` |
| [Li & Horvitz 1953](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1716461/)| moments-based | `LiHorvitz` |
| [Loiselle 1995](https://bsapubs.onlinelibrary.wiley.com/doi/abs/10.1002/j.1537-2197.1995.tb12679.x) | moments-based | `Loiselle` |
| [Lynch 1988](https://pubmed.ncbi.nlm.nih.gov/3193879/) | moments-based | `Lynch` |
| [Lynch/Li 1993](https://pubmed.ncbi.nlm.nih.gov/8514326/) | moments-based | `LynchLi` |
| [Lynch & Ritland 1999](https://www.genetics.org/content/152/4/1753.short) | moments-based | `LynchRitland` |
| [Moran 1950](https://www.jstor.org/stable/2332142?origin=crossref&seq=1#metadata_info_tab_contents) | moments-based | `Moran` |
| [Queller & Goodnight 1989](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1989.tb04226.x) | moments-based | `QuellerGoodnight` |
| [Ritland 1996](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1558-5646.1996.tb02347.x) | moments-based | `Ritland` |
| [Milligan 2003](https://pubmed.ncbi.nlm.nih.gov/12663552/) "DyadML" | maximum-likelihood | incomplete* |
| [Wang 2002](https://www.genetics.org/content/160/3/1203.short) | moments-based | incomplete* |

:::note *more kinship estimators
Contact us or submit a pull request if you're interested in contributing to the kinship methods. DyadML and Wang (2002) estimators are currently incomplete and the others
could use some optimizations. More help is always welcomed! Our wishlist also includes the KING method :smile:
:::

## Posthoc analyses
There are several different kinds of things you can do with kinship information (e.g. network analysis), and one that's provided is lovingly
called `kinshipposthoc()`, which performs a permutation analysis to
test if within-population relatedness is significantly greater than between-population relatedness for each population. We recommend that you
correct for multiple testing using `MultipleTesting.jl`.

```julia
kinshipposthoc(data::PopData, results::Union{DataFrame, NamedTuple}; iterations::Int)
```
#### Arguments
- `data` : A PopData object
- `results` : the DataFrame or NamedTuple results from `kinship()`

#### Keyword Arguments
- `iterations` : number of iterations for the permutation tests (default: `20000`)

:::tip not a great name
We admit "kinshipposthoc" is not a great name for this function. Please
contact us with better ideas! :grin:
:::

**Example**
```
julia> cats = @nancycats ;

julia> rel_out = kinship(cats, method = [Ritland, Moran], iterations = 100);

julia> kinshipposthoc(cats, rel_out)
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
The kinship methods were dutifully written and verified against R analogues by Jason Selwyn. These anaylses can take a while, especially if performing bootstrapping; we provide a progress bar via [Term.jl](https://github.com/FedeClaudi/Term.jl) so you can move on and focus on other things in the meantime. 
