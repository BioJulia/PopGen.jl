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
pairwisefst(data::PopData; method::String, iterations::Int64)
```
Calculate pairwise $F_{ST}$ between populations in a `PopData` object. Set
`iterations` to a value greater than `0` to perform a one-tailed permutation
test to obtain p-values of statistical significance. The permutations respect the
original population sizes. The pairwise $F_{ST}$ values will appear below the diagonal
 (which all be zeros), and the corresponding p-values will appear above the diagonal. 
 A progress bar is provided for significance testing.

:::note custom output type
The returned object is a custom `PairwiseFST` type with the fields `results` (stores the dataframe of $F_{ST}$ values) and `method` (a string of which method was used to calculate it). This was done to define a custom `show` method to make the results a little nicer, and so you never lose track of which method was used for the calculation. If you want to access the dataframe directly, you will need to do so with `varname.results` where `varname` is whatever you named the output.  
:::

### Arguments
- `data::PopData`: a PopData object you wish to perform the calculation on

### Keyword Arguments
- `method::String`: which $F_{ST}$ calculation method you would like to use
    - `"Hudson92"`: the [Hudson et al. (1992)](https://www.genetics.org/content/132/2/583) method (only for biallelic data)
    - `"WC84"`: the [Weir & Cockerham (1984)](https://www.jstor.org/stable/2408641?casa_token=_0gGbCbYpqMAAAAA:f9BvW9Xvx_8WaWSaRN4iqg0HB7KkaP21712ds28cTjhsvVQrYRTyHon7hKCcyHLcmTRA9H_1oM5iF3TZAl5xPm5gil2GmcGzHyEFFYAOl8pDVEBMQQ&seq=1#metadata_info_tab_contents) method (default)
    - `"Nei87"`: [Nei's (1987)](https://books.google.com/books?hl=en&lr=&id=UhRSsLkjxDgC&oi=fnd&pg=PP11&ots=Qu7vO8EMmw&sig=T6cTISYEEm-hL8aWU8EgeGgzP5E#v=onepage&q&f=false) genetic distance method
- `iterations::Int64`: the number of iterations for signficance testing (default: `0`)

**Examples**

<Tabs
  block={true}
  defaultValue="wop"
  values={[
    { label: 'Without significance testing', value: 'wop', },
    { label: 'With significance testing', value: 'wp', },
  ]
}>
<TabItem value="wop">

```
julia> sharks = @gulfsharks ;

julia> pairwisefst(sharks, method = "WC84")
Pairwise FST: Weir & Cockerham
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
julia> cats = @nancycats ;

julia> pairwisefst(cats, iterations = 500)

Below diagonal: FST values | Above diagonal: P values
Pairwise FST: Weir & Cockerham (with p-values)
    │ 1        2        3        4        5        6        7        8        9        10       11       12       13       14       15       16       17    
────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
  1 │ 0.0      0.002    0.002    0.002    0.002    0.002    0.002    0.004    0.002    0.004    0.008    0.002    0.002    0.002    0.002    0.002    0.006
  2 │ 0.13077  0.0      0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002
  3 │ 0.0846   0.13112  0.0      0.03     0.002    0.03     0.002    0.002    0.002    0.002    0.014    0.002    0.002    0.002    0.002    0.002    0.002
  4 │ 0.0714   0.10729  0.01939  0.0      0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002
  5 │ 0.06928  0.12973  0.06186  0.04689  0.0      0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002
  6 │ 0.09329  0.09765  0.02451  0.03816  0.07535  0.0      0.002    0.002    0.002    0.002    0.02     0.002    0.002    0.002    0.002    0.002    0.002
  7 │ 0.07283  0.10678  0.10116  0.07909  0.11943  0.0764   0.0      0.008    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002
  8 │ 0.05247  0.08163  0.06064  0.03795  0.08567  0.05967  0.04643  0.0      0.002    0.004    0.004    0.002    0.002    0.002    0.002    0.002    0.002
  9 │ 0.11191  0.12121  0.10559  0.073    0.1279   0.09913  0.14718  0.11268  0.0      0.002    0.002    0.002    0.002    0.002    0.002    0.002    0.002
 10 │ 0.05744  0.07593  0.09633  0.0975   0.09325  0.08429  0.10651  0.05301  0.1033   0.0      0.008    0.002    0.002    0.002    0.004    0.002    0.002
 11 │ 0.05164  0.09153  0.02871  0.02973  0.05785  0.0288   0.08206  0.05314  0.05899  0.04586  0.0      0.002    0.002    0.002    0.002    0.002    0.002
 12 │ 0.07807  0.11885  0.08425  0.05782  0.09224  0.08536  0.07691  0.07919  0.0865   0.11416  0.07992  0.0      0.002    0.002    0.002    0.002    0.002
 13 │ 0.1012   0.12755  0.11318  0.09985  0.11603  0.09097  0.13132  0.1041   0.15204  0.12559  0.10167  0.11873  0.0      0.002    0.002    0.002    0.002
 14 │ 0.07518  0.10531  0.04524  0.04166  0.07639  0.06064  0.10501  0.04941  0.07358  0.0866   0.04558  0.07863  0.09973  0.0      0.002    0.002    0.002
 15 │ 0.06767  0.09057  0.06297  0.05262  0.09979  0.08144  0.10854  0.0442   0.07541  0.05113  0.05764  0.10327  0.1044   0.0558   0.0      0.002    0.002
 16 │ 0.11927  0.1119   0.1034   0.0772   0.09671  0.0559   0.12758  0.10547  0.09288  0.1077   0.06333  0.11881  0.12398  0.0841   0.0972   0.0      0.002
 17 │ 0.06563  0.14166  0.10072  0.08063  0.07786  0.12416  0.09867  0.11597  0.12334  0.10998  0.06634  0.06166  0.12915  0.08636  0.09778  0.15081  0.0
```

</TabItem>
</Tabs>