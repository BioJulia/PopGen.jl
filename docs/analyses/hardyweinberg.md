---
id: hardyweinberg
title: Hardy-Weinberg Equilibrium
sidebar_label: Hardy-Weinberg Equilibrium
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />

Testing for Hardy-Weinberg Equilibrium (often abbreviated _HWE_) is a fairly common practice in population genetics. In a two-allele system, the HWE equation is defined as: $p^2 + 2pq + q^2 = 1$ , where $p$ is the frequency of the first allele and $q$ is the frequency of the second allele. The formula describes the frequency of all possible genotypes where

| HWE variable | Genotype |    State     |
| :----------: | :------: | :----------: |
|    $p^2$    |   "pp"   |  homozygous  |
|    $q^2$    |   "qq"   |  homozygous  |
|    $2pq$     |   "pq"   | heterozygous |

Testing for deviation from HWE is usually done with a Chi-Squared ($\chi^2$) test, where one compares the observed genotype frequencies to the expected genotype frequencies given the observed allele frequencies at a locus. Specifically, the equation is
$$
\sum{\frac{(observed - expected)^2}{expected}}
$$
where $observed$ is the observed genotype frequency and $expected$ is the expected genotype frequency for a locus. To generate our test statistic, we calculate the degrees of freedom: 
$$
degrees\ of\ freedom = n_{expected\ genotypes} - n_{observed\ alleles}
$$ 
and use this as the parameter for our $\chi^2$ distribution, followed by a cumulative density function using this $\chi^2$ distribution and our $\chi^2$ value calculated above.

## Chi-Squared Test

```julia
hwe_test(data::PopData, by::String = "locus", correction::String = "none")
```

Calculate chi-squared test of HWE for each locus and returns observed and expected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use `by = "population"` to perform this separately for each population (default: `"locus"`). Use `correction =` to specify a P-value correction method for multiple testing (recommended). For convenience, the correction method is appended to the name of the column, so you will always know how those P-values were adjusted.

### arguments

- `data` : the input `PopData`
- `by_pop =` : `false` (default) or `true` for hwe-by-population
- `correction =`  : a string specifying a P-value adjustment type (default: "none")

### `correction` methods

- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"`  : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment

:thinking: For more information on multiple testing adjustments, see [MultipleTesting.jl](https://juliangehring.github.io/MultipleTesting.jl/stable/)

### example
<Tabs
  block={true}
  defaultValue="hwe"
  values={[
    { label: 'HWE Chi-Sq', value: 'hwe', },
    { label: 'HWE with P adjustment', value: 'hwe_p', },
    { label: 'HWE with by population', value: 'hwe_pop', },
  ]
}>
<TabItem value="hwe">

```
julia> hwe_test(gulfsharks())
2213×4 DataFrame
│ Row  │ locus        │ chisq     │ df    │ P        │
│      │ Categorical… │ Float64   │ Int64 │ Float64? │
├──────┼──────────────┼───────────┼───────┼──────────┤
│ 1    │ contig_35208 │ 0.690979  │ 3     │ 0.875324 │
│ 2    │ contig_23109 │ 1.71251   │ 1     │ 0.19066  │
│ 3    │ contig_4493  │ 2.0785    │ 1     │ 0.149387 │
│ 4    │ contig_10742 │ 0.249703  │ 1     │ 0.617284 │
│ 5    │ contig_14898 │ 2.10859   │ 1     │ 0.146474 │
│ 6    │ contig_8483  │ 0.019229  │ 1     │ 0.889712 │
⋮
│ 2207 │ contig_18959 │ 0.968083  │ 1     │ 0.325159 │
│ 2208 │ contig_43517 │ 0.655517  │ 1     │ 0.418147 │
│ 2209 │ contig_27356 │ 1.71412   │ 1     │ 0.190451 │
│ 2210 │ contig_475   │ 0.0754715 │ 1     │ 0.783531 │
│ 2211 │ contig_19384 │ 2.10124   │ 1     │ 0.147179 │
│ 2212 │ contig_22368 │ 0.473923  │ 1     │ 0.491188 │
│ 2213 │ contig_2784  │ 0.0452162 │ 1     │ 0.831607 │
```

</TabItem>
<TabItem value="hwe_p">

```
julia> hwe_test(gulfsharks(), correction = "bh")
2213×5 DataFrame
│ Row  │ locus        │ chisq     │ df    │ P        │ P_bh     │
│      │ Categorical… │ Float64   │ Int64 │ Float64? │ Float64? │
├──────┼──────────────┼───────────┼───────┼──────────┼──────────┤
│ 1    │ contig_35208 │ 0.690979  │ 3     │ 0.875324 │ 0.999911 │
│ 2    │ contig_23109 │ 1.71251   │ 1     │ 0.19066  │ 0.999911 │
│ 3    │ contig_4493  │ 2.0785    │ 1     │ 0.149387 │ 0.999911 │
│ 4    │ contig_10742 │ 0.249703  │ 1     │ 0.617284 │ 0.999911 │
│ 5    │ contig_14898 │ 2.10859   │ 1     │ 0.146474 │ 0.999911 │
│ 6    │ contig_8483  │ 0.019229  │ 1     │ 0.889712 │ 0.999911 │
⋮
│ 2207 │ contig_18959 │ 0.968083  │ 1     │ 0.325159 │ 0.999911 │
│ 2208 │ contig_43517 │ 0.655517  │ 1     │ 0.418147 │ 0.999911 │
│ 2209 │ contig_27356 │ 1.71412   │ 1     │ 0.190451 │ 0.999911 │
│ 2210 │ contig_475   │ 0.0754715 │ 1     │ 0.783531 │ 0.999911 │
│ 2211 │ contig_19384 │ 2.10124   │ 1     │ 0.147179 │ 0.999911 │
│ 2212 │ contig_22368 │ 0.473923  │ 1     │ 0.491188 │ 0.999911 │
│ 2213 │ contig_2784  │ 0.0452162 │ 1     │ 0.831607 │ 0.999911 │
```

</TabItem>
<TabItem value="hwe_pop">

```
julia> hwe_test(gulfsharks(), by = "population")
15491×5 DataFrame
│ Row   │ locus        │ population     │ chisq       │ df    │ P        │
│       │ Categorical… │ Categorical…   │ Float64     │ Int64 │ Float64? │
├───────┼──────────────┼────────────────┼─────────────┼───────┼──────────┤
│ 1     │ contig_35208 │ Cape Canaveral │ 1.01787     │ 1     │ 0.313025 │
│ 2     │ contig_35208 │ Georgia        │ 0.725926    │ 1     │ 0.394207 │
│ 3     │ contig_35208 │ South Carolina │ 0.571429    │ 1     │ 0.449692 │
│ 4     │ contig_35208 │ Florida Keys   │ 0.276208    │ 3     │ 0.964439 │
│ 5     │ contig_35208 │ Mideast Gulf   │ 0.598475    │ 1     │ 0.439161 │
│ 6     │ contig_35208 │ Northeast Gulf │ 0.519382    │ 3     │ 0.914613 │
⋮
│ 15485 │ contig_2784  │ Cape Canaveral │ 2.9277e-11  │ 0     │ missing  │
│ 15486 │ contig_2784  │ Georgia        │ 1.05242e-11 │ 0     │ missing  │
│ 15487 │ contig_2784  │ South Carolina │ 1.94026e-11 │ 0     │ missing  │
│ 15488 │ contig_2784  │ Florida Keys   │ 0.0374777   │ 1     │ 0.846496 │
│ 15489 │ contig_2784  │ Mideast Gulf   │ 0.0399408   │ 1     │ 0.841596 │
│ 15490 │ contig_2784  │ Northeast Gulf │ 1.47338e-11 │ 0     │ missing  │
│ 15491 │ contig_2784  │ Southeast Gulf │ 0.0138787   │ 1     │ 0.90622  │
```
When doing this test by population, you may notice some loci have `missing` P-values for certain populations, indicating that this locus is missing for that population. 

</TabItem>
</Tabs>

## Interpreting the results
Since the results are in table form, you can easily process the table using `DataFramesMeta.jl` or `Query.jl` to find loci above or below the alpha threshold you want. As an example, let's perform an HWE-test on the `nancycats` data without any P-value adjustments:
```julia
julia> ncats_hwe = hwe_test(nancycats() , by = "population") ;
```
Now, we can use `DataFramesMeta.jl` to easily filter this table and leave only what we're interested in:
```julia
using DataFramesMeta

julia> @where(ncats_hwe, :P .!== missing, :P .<= 0.05)
```
With this command, we invoke the `@where` filtering macro, then specify our `ncats_hwe` table, the `:P` column of P-values, and then specify two filtering parameters: 
1. the P-values are not `missing`
2. the P-values are less than or equal to 0.05. 

Doing this results in a table that now only includes non-missing P-values of 0.05 or lower:
```
46×5 DataFrame
│ Row │ locus │ population │ chisq   │ df    │ P           │
│     │ Cat…  │ Cat…       │ Float64 │ Int64 │ Float64?    │
├─────┼───────┼────────────┼─────────┼───────┼─────────────┤
│ 1   │ fca8  │ 2          │ 29.8624 │ 15    │ 0.0124275   │
│ 2   │ fca8  │ 3          │ 39.72   │ 21    │ 0.00804172  │
│ 3   │ fca8  │ 4          │ 64.3039 │ 45    │ 0.0308567   │
│ 4   │ fca8  │ 6          │ 34.0525 │ 21    │ 0.0357734   │
│ 5   │ fca8  │ 7          │ 37.6531 │ 15    │ 0.00101511  │
│ 6   │ fca23 │ 1          │ 26.0494 │ 10    │ 0.00367434  │
⋮
│ 40  │ fca96 │ 14         │ 58.0854 │ 36    │ 0.0112879   │
│ 41  │ fca96 │ 15         │ 26.51   │ 15    │ 0.0329918   │
│ 42  │ fca96 │ 2          │ 92.0716 │ 15    │ 4.06675e-13 │
│ 43  │ fca37 │ 1          │ 10.1562 │ 3     │ 0.0172836   │
│ 44  │ fca37 │ 3          │ 13.4815 │ 6     │ 0.0359962   │
│ 45  │ fca37 │ 5          │ 30.0206 │ 6     │ 3.89561e-5  │
│ 46  │ fca37 │ 8          │ 33.5111 │ 21    │ 0.0408481   │
```

### Visualizing the results
While not strictly necessary, it might sometimes make sense to generate of heatmap of the results for easier visualization. This is feasible for the `nancycats` data, but when loci are in the hundreds or thousands, this method quickly becomes counterproductive. In any case, here is a simple example of the HWE results for `nancycats` plotted as a heatmap using [VegaLite.jl](https://github.com/queryverse/VegaLite.jl):
```julia
using VegaLite

julia> ncats_hwe = hwe_test(nancycats() , by = "population", correction = "bonferroni");

julia> ncats |> @vlplot(:rect, :locus, :population, color=:P_bonferroni)
```
![hwe_test](/PopGen.jl/img/hwe_test.png)

