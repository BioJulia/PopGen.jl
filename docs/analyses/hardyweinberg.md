# Hardy-Weinberg Equilibrium

Testing for Hardy-Weinberg Equilibrium (often abbreviated _HWE_) is a fairly common practice in population genetics. In a two-allele system, the HWE equation is defined as:
$$p^2 + 2pq + q^2 = 1$$
where $p$ is the frequency of the first allele and $q$ is the frequency of the second allele. The formula describes the frequency of all possible genotypes where

| HWE variable | Genotype |    State     |
| :----------: | :------: | :----------: |
|    $p^2$    |   "pp"   |  homozygous  |
|    $q^2$    |   "qq"   |  homozygous  |
|    $2pq$     |   "pq"   | heterozygous |

Testing for deviation from HWE is usually done with a Chi-Squared test, where one compares the observed genotype frequencies to the expected genotype frequencies given the observed allele frequencies at a locus. Specifically the equation is
$$\sum{\frac{(observed - expected)^2}{expected}}$$
where $observed$ is the observed genotype frequency and $expected$ is the expected genotype frequency for a locus. To generate our test statistic, we calculate the degrees of freedom: 
$$degrees\ of\ freedom = n_{expected\ genotypes} - n_{observed\ alleles}$$ 
and use this as the parameter for our Chi Squared distribution, followed by a cumulative density function using this Chi Squared distribution and our Chi-Squared value calculated above.
## Chi-Squared Test

```julia
hwe_test(x::PopObj, by_pop::Bool = false, correction::String = "none")
```

Calculate chi-squared test of HWE for each locus and returns observed and expected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use `by_pop = true` to perform this separately for each population (default: by_pop = false) and return a NamedTuple with the names corresponding to the population names. Use `correction =` to specify a P-value correction method for multiple testing (recommended).

### arguments

- `x` : the input `PopObj`
- `by_pop =` : `false` (default) or `true` for hwe-by-population
- `correction =`  : a string specifying a P-value adjustment type (default: "none")

### `correction` methods

- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` or `"b-h"` : Benjamini-Hochberg adjustment
- `"by"` or `"b-y"`: Benjamini-Yekutieli adjustment
- `"bl"` or `"b-l"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forward stop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` or `"b-c"` : Barber-Candès adjustment

:thinking: For more information on multiple testing adjustments, see [MultipleTesting.jl](https://juliangehring.github.io/MultipleTesting.jl/stable/)

### example
:::: tabs card true
::: tab HWE Chi-Squared example
```julia
hwe_test(gulfsharks(), correction = "bh")
```
:::
::: tab output
```
2213×7 DataFrames.DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²        │ DF       │ P          │ Pcorr      │
│      │ String       │ Float64   │ Float64   │ Float64⍰  │ Float64⍰ │ Float64⍰   │ Float64⍰   │
├──────┼──────────────┼───────────┼───────────┼───────────┼──────────┼────────────┼────────────┤
│ 1    │ contig_35208 │ 0.419811  │ 0.398051  │ 0.690981  │ 3.0      │ 0.875323   │ 0.948912   │
│ 2    │ contig_23109 │ 0.262136  │ 0.288434  │ 1.71251   │ 1.0      │ 0.19066    │ 0.295543   │
│ 3    │ contig_4493  │ 0.205742  │ 0.228532  │ 2.0785    │ 1.0      │ 0.149387   │ 0.234857   │
│ 4    │ contig_10742 │ 0.0666667 │ 0.0644444 │ 0.249703  │ 1.0      │ 0.617284   │ 0.807471   │
│ 5    │ contig_14898 │ 0.240566  │ 0.21875   │ 2.10859   │ 1.0      │ 0.146475   │ 0.230442   │
│ 6    │ contig_8483  │ 0.0188679 │ 0.0186899 │ 0.019229  │ 1.0      │ 0.889712   │ 0.948912   │
│ 7    │ contig_8065  │ 0.0801887 │ 0.0769736 │ 0.369866  │ 1.0      │ 0.543077   │ 0.732045   │
│ 8    │ contig_14708 │ 0.0616114 │ 0.0597134 │ 0.213168  │ 1.0      │ 0.644295   │ 0.828219   │
│ 9    │ contig_2307  │ 0.0289855 │ 0.0285654 │ 0.0447664 │ 1.0      │ 0.832434   │ 0.938862   │
│ 10   │ contig_14564 │ 0.2       │ 0.209751  │ 0.453809  │ 1.0      │ 0.500532   │ 0.6894     │
│ 11   │ contig_15269 │ 0.15566   │ 0.146505  │ 1.51012   │ 3.0      │ 0.679938   │ 0.857489   │
│ 12   │ contig_24796 │ 0.161137  │ 0.201613  │ 8.50406   │ 1.0      │ 0.00354355 │ 0.00596524 │
│ 13   │ contig_14251 │ 0.490566  │ 0.49782   │ 0.0450074 │ 1.0      │ 0.83199    │ 0.938862   │
│ 14   │ contig_44797 │ 0.0240385 │ 0.0237495 │ 0.0307836 │ 1.0      │ 0.860724   │ 0.944468   │
│ 15   │ contig_43681 │ 0.42381   │ 0.472778  │ 2.25286   │ 1.0      │ 0.133368   │ 0.211628   │
│ 16   │ contig_24115 │ 0.0333333 │ 0.0327778 │ 0.0603275 │ 1.0      │ 0.805979   │ 0.931567   │
│ 17   │ contig_5456  │ 0.0471698 │ 0.0460573 │ 0.12369   │ 1.0      │ 0.725066   │ 0.892017   │
│ 18   │ contig_21698 │ 0.0758294 │ 0.0729543 │ 0.327695  │ 1.0      │ 0.567019   │ 0.758715   │
⋮
│ 2195 │ contig_8479  │ 0.0235849 │ 0.0598856 │ 212.031   │ 3.0      │ 0.0        │ 0.0        │
│ 2196 │ contig_47462 │ 0.45283   │ 0.502803  │ 213.589   │ 3.0      │ 0.0        │ 0.0        │
│ 2197 │ contig_4095  │ 0.259434  │ 0.275398  │ 212.0     │ 3.0      │ 0.0        │ 0.0        │
│ 2198 │ contig_7239  │ 0.0660377 │ 0.0729352 │ 212.25    │ 3.0      │ 0.0        │ 0.0        │
│ 2199 │ contig_40507 │ 0.0330189 │ 0.0688746 │ 212.062   │ 3.0      │ 0.0        │ 0.0        │
│ 2200 │ contig_42145 │ 0.358491  │ 0.31754   │ 5.95351   │ 3.0      │ 0.113894   │ 0.182428   │
│ 2201 │ contig_1033  │ 0.363208  │ 0.352872  │ 0.181856  │ 1.0      │ 0.669783   │ 0.85219    │
│ 2202 │ contig_2798  │ 0.122642  │ 0.140219  │ 212.952   │ 3.0      │ 0.0        │ 0.0        │
│ 2203 │ contig_12991 │ 0.0518868 │ 0.112685  │ 212.161   │ 3.0      │ 0.0        │ 0.0        │
│ 2204 │ contig_22981 │ 0.188679  │ 0.200783  │ 0.770434  │ 1.0      │ 0.380083   │ 0.551307   │
│ 2205 │ contig_15342 │ 0.268868  │ 0.29084   │ 212.328   │ 3.0      │ 0.0        │ 0.0        │
│ 2206 │ contig_24711 │ 0.273585  │ 0.288136  │ 0.540694  │ 1.0      │ 0.462145   │ 0.645777   │
│ 2207 │ contig_18959 │ 0.466981  │ 0.437422  │ 0.968086  │ 1.0      │ 0.325158   │ 0.482394   │
│ 2208 │ contig_43517 │ 0.103774  │ 0.150454  │ 212.675   │ 3.0      │ 0.0        │ 0.0        │
│ 2209 │ contig_27356 │ 0.0660377 │ 0.0906016 │ 213.73    │ 3.0      │ 0.0        │ 0.0        │
│ 2210 │ contig_475   │ 0.367925  │ 0.375     │ 0.0754717 │ 1.0      │ 0.78353    │ 0.926932   │
│ 2211 │ contig_19384 │ 0.0613208 │ 0.11264   │ 214.152   │ 3.0      │ 0.0        │ 0.0        │
│ 2212 │ contig_22368 │ 0.0896226 │ 0.11224   │ 212.481   │ 3.0      │ 0.0        │ 0.0        │
│ 2213 │ contig_2784  │ 0.0283019 │ 0.0908241 │ 212.047   │ 3.0      │ 0.0        │ 0.0        │
```
:::
::: tab HWE by locus per population
```
julia> hwe_test(gulfsharks(), by_pop = true, correction = "bh")

```
:::
::::