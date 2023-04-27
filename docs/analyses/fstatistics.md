---
id: fstatistics
title: Pairwise F-Statistics
sidebar_label: Pairwise F-Statistics
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />


## Background
Originating with [Wright's $F$ fixation index](https://www.genetics.org/content/16/2/97) (inbreeding coefficient), $F$ has 
since evolved into a whole slew of statistics used to describe the differentiation/divergence within or between groups. As
you may have seen from `summary()`, there are a common group of these indices ($F_{IS}$, $F_{IT}$, and $F_{ST}$) that compare the $F$ at various hierarchical levels. The notation is pretty straightforward: $I$ is for individuals, $T$ is the total, and $S$ is for subpopulations.

| F-statistic  | Compares           |   Against         |
| :----------: |     :------       |   :----------    |
|    $F_{IS}$  |   **I**ndividual   | **S**ubpopulation |
|    $F_{IT}$  |   **I**ndividual   | **T**otal         |
|    $F_{ST}$  |  **S**ubpopulation | **T**otal         |


Often, we are interested in pairwise $F_{ST}$, which is a type of coefficient that helps us infer how panmictic (fully mixed) two groups of interest are. A colloquial way of phrasing that is "how much genetic mixing is there between these two groups?". The value of $F_{ST}$ (and its derivatives) should typically range between 0 and 1 and can be interpreted like so:

| $F_{ST}$ value    |                  Interpretation           |
| :----------:      |           :-----------------------:       |
|         0         |   the two groups are completely panmictic |
|         1         |   the two groups are completely isolated  |

However, it's not a linear relationship, and Wright considered 0.125 as the cutoff for when to determine populations as divergent.

:::note $F_{ST}$ isn't everything
An important caveat to always consider is that $F_{ST}$ is just one tool to help us understand trends and not the entire picture.
The genetic data we collect is just a snapshot in current time and populations can be completely isolated but still have near-zero
$F_{ST}$ values for a number of reasons (slow divergence time, recent introgression, etc.). Significance testing helps add context
to observed $F_{ST}$ values.
:::

## Pairwise FST
```julia
pairwisefst(data::PopData; method::Function, by::String = "global", iterations::Int64)
```
Calculate pairwise $F_{ST}$ between populations in a `PopData` object. Set `iterations` 
to a value greater than `0` to perform a single-tailed permutation test to obtain
P-values of statistical significance. Use `by = "locus"` to perform a locus-by-locus FST for
population pairs (iterations and significance testing ignored). `WeirCockerham is not yet implmented for by-locus $F_{ST}$. 

:::note custom output type
The returned object for is a custom `PairwiseFST` type with the fields `results` (stores the dataframe of $F_{ST}$ values) and `method` (a string of which method was used to calculate it). This was done to define a custom `show` method to make the results a little nicer, and so you never lose track of which method was used for the calculation. If you want to access the dataframe directly, you will need to do so with `varname.results` where `varname` is whatever you named the output.  
:::

### Arguments
- `data::PopData`: a PopData object you wish to perform the calculation on

### Keyword Arguments
- `method::Function`: which $F_{ST}$ calculation method you would like to use
    - `AMOVA`: the Analysis of Molecular Variance method from [Bird et al. 2011](https://www.researchgate.net/publication/229089010_Detecting_and_measuring_genetic_differentiation) (default)
    - `Hudson`: the [Hudson et al. (1992)](https://www.genetics.org/content/132/2/583) method (only for biallelic data)
    - `WeirCockerham`: the [Weir & Cockerham (1984)](https://www.jstor.org/stable/2408641?casa_token=_0gGbCbYpqMAAAAA:f9BvW9Xvx_8WaWSaRN4iqg0HB7KkaP21712ds28cTjhsvVQrYRTyHon7hKCcyHLcmTRA9H_1oM5iF3TZAl5xPm5gil2GmcGzHyEFFYAOl8pDVEBMQQ&seq=1#metadata_info_tab_contents) method
    - `Nei`: [Nei (1987)](https://books.google.com/books?hl=en&lr=&id=UhRSsLkjxDgC&oi=fnd&pg=PP11&ots=Qu7vO8EMmw&sig=T6cTISYEEm-hL8aWU8EgeGgzP5E#v=onepage&q&f=false) genetic distance method
- `by::String`: perfrom a `"global"` pairwise $F_{ST}$ or `"locus"` for locus-by-locus (ignores significance testing)
- `iterations::Int64`: the number of iterations for signficance testing (default: `0`)

**Examples**

<Tabs
  block={true}
  defaultValue="wop"
  values={[
    { label: 'without significance testing', value: 'wop', },
    { label: 'with significance testing', value: 'wp', },
    { label: 'by locus', value: 'loc', },
  ]
}>
<TabItem value="wop">

```julia
julia> sharks = @gulfsharks ;

julia> pairwisefst(sharks, method = WeirCockerham)
Pairwise FST: WeirCockerham
                CapeCanaveral  Georgia   SouthCarolina  FloridaKeys  MideastGulf  NortheastGulf  SoutheastGulf 
───────────────────────────────────────────────────────────────────────────────────────────────────────────────
 CapeCanaveral        0.0       0.0            0.0          0.0          0.0            0.0                0.0
       Georgia        0.00081   0.0            0.0          0.0          0.0            0.0                0.0
 SouthCarolina       -0.0003   -0.00076        0.0          0.0          0.0            0.0                0.0
   FloridaKeys        0.00282   0.00202        0.00204      0.0          0.0            0.0                0.0
   MideastGulf        0.00423   0.00354        0.00329      0.00042      0.0            0.0                0.0
 NortheastGulf        0.00264   0.00147        0.00146     -7.0e-5      -0.00023        0.0                0.0
 SoutheastGulf        0.00312   0.00222        0.00191     -3.0e-5       0.00079        0.00118            0.0
```

</TabItem>
<TabItem value="wp">

```julia
julia> sharks = @gulfsharks ;

julia> pairwisefst(sharks, iterations = 100)

Below diagonal: FST values | Above diagonal: P values
Pairwise FST: WeirCockerham (with p-values)
                CapeCanaveral  Georgia   SouthCarolina  FloridaKeys  MideastGulf  NortheastGulf  SoutheastGulf 
───────────────────────────────────────────────────────────────────────────────────────────────────────────────
 CapeCanaveral        0.0       0.17           0.66         0.01         0.01           0.08              0.01
       Georgia        0.00081   0.0            0.83         0.01         0.01           0.05              0.02
 SouthCarolina       -0.0003   -0.00076        0.0          0.01         0.01           0.09              0.03
   FloridaKeys        0.00282   0.00202        0.00204      0.0          0.23           0.6               0.65
   MideastGulf        0.00423   0.00354        0.00329      0.00042      0.0            0.64              0.18
 NortheastGulf        0.00264   0.00147        0.00146     -7.0e-5      -0.00023        0.0               0.15
 SoutheastGulf        0.00312   0.00222        0.00191     -3.0e-5       0.00079        0.00118           0.0
```

</TabItem>
<TabItem value="loc">

```julia
julia> sharks = @gulfsharks ;
julia> pairwisefst(sharks, method = Nei, by = "locus")
Pairwise FST: Nei
pop1           pop2           locus         fst      
String         String         String        Float64  
─────────────────────────────────────────────────────
CapeCanaveral  Georgia        contig_35208  -0.01309
CapeCanaveral  Georgia        contig_23109  -0.0235
CapeCanaveral  Georgia        contig_4493   -0.02533
CapeCanaveral  Georgia        contig_10742   0.01365
CapeCanaveral  Georgia        contig_14898   0.00107
CapeCanaveral  Georgia        contig_8483   -0.00372
      ⋮              ⋮             ⋮           ⋮
NortheastGulf  SoutheastGulf  contig_43517  -0.00146
NortheastGulf  SoutheastGulf  contig_27356   0.02708
NortheastGulf  SoutheastGulf  contig_475     0.0081
NortheastGulf  SoutheastGulf  contig_19384   0.05055
NortheastGulf  SoutheastGulf  contig_22368  -0.0036
NortheastGulf  SoutheastGulf  contig_2784   -0.00066
                                            46377 rows omitted
```

</TabItem>
</Tabs>