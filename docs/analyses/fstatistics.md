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
pairwise_fst(data::PopData; method::String)
```
Calculate pairwise $F_{ST}$ between populations in a `PopData` object.

:::note custom output type
The returned object is a custom `PairwiseFST` type with the fields `results` (stores the dataframe of $F_{ST}$ values) and `method` (a string of which method was used to calculate it). This was done to define a custom `show` method to make the results a little nicer, and so you never lose track of which method was used for the calculation. If you want to access the dataframe directly, you will need to do so with `varname.results`.  
:::

### Arguments
- `data::PopData`: a PopData object you wish to perform the calculation on

### Keyword Arguments
- `method::String`: which $F_{ST}$ calculation method you would like to use
    - `"WC84"`: the [Weir & Cockerham (1984)](https://www.jstor.org/stable/2408641?casa_token=_0gGbCbYpqMAAAAA:f9BvW9Xvx_8WaWSaRN4iqg0HB7KkaP21712ds28cTjhsvVQrYRTyHon7hKCcyHLcmTRA9H_1oM5iF3TZAl5xPm5gil2GmcGzHyEFFYAOl8pDVEBMQQ&seq=1#metadata_info_tab_contents) method (default)
    - `"Nei87"`: [Nei's (1987)](https://books.google.com/books?hl=en&lr=&id=UhRSsLkjxDgC&oi=fnd&pg=PP11&ots=Qu7vO8EMmw&sig=T6cTISYEEm-hL8aWU8EgeGgzP5E#v=onepage&q&f=false) genetic distance method

**Example**
```
julia> sharks = @gulfsharks ;

julia> pairwise_fst(sharks, method = "WC84")
Pairwise FST: Weir & Cockerham
                │ Cape Canaveral  Georgia   South Carolina  Florida Keys  Mideast Gulf  Northeast Gulf  Southeast Gulf 
────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────
 Cape Canaveral │        0.0       0.0             0.0           0.0           0.0             0.0                 0.0
        Georgia │        0.00081   0.0             0.0           0.0           0.0             0.0                 0.0
 South Carolina │       -0.0003   -0.00076         0.0           0.0           0.0             0.0                 0.0
   Florida Keys │        0.00282   0.00202         0.00204       0.0           0.0             0.0                 0.0
   Mideast Gulf │        0.00423   0.00354         0.00329       0.00042       0.0             0.0                 0.0
 Northeast Gulf │        0.00264   0.00147         0.00146      -7.0e-5       -0.00023         0.0                 0.0
 Southeast Gulf │        0.00312   0.00222         0.00191      -3.0e-5        0.00079         0.00118             0.0

```