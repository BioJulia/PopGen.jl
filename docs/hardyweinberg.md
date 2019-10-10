Testing for Hardy-Weinberg Equilibrium (often abbreviated to "HW" or "HWE") is a fairly common practice in population genetics. In a two-allele system (alleles *p* and *q*), the HWE equation is defined as *p^2^ + 2pq + q^2^ = 1*. 

Testing for deviation from HWE is usually done with a Chi-Squared test, where one compares the observed genotype frequencies to the expected genotype frequencies given the observed allele frequencies at a locus. 

## Test for Hardy-Weinberg equilibrium

```julia
hwe_test(x::PopObj, by_pop::Bool = false correction::{String} = "none")
```

Calculate chi-squared test of HWE for each locus and returns observed and expected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use `by_pop = true` to perform this separately for each population (default: by_pop = false) and return a NamedTuple with the names corresponding to the population names. Use `correction =` to specify a P-value
correction method for multiple testing.

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

```julia tab="example"
hwe_test(gulfsharks(), correction = "bh")
```

``` tab="output"
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

``` tab="by_pop = true"
julia> hwe_test(gulfsharks(), by_pop = true, correction = "bh")

(Cape Canaveral = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²        │ DF       │ P           │ Pcorr       │
│      │ String       │ Float64⍰  │ Float64⍰  │ Float64⍰  │ Float64⍰ │ Float64⍰    │ Float64⍰    │
├──────┼──────────────┼───────────┼───────────┼───────────┼──────────┼─────────────┼─────────────┤
│ 1    │ contig_35208 │ 0.333333  │ 0.427438  │ 1.01787   │ 1.0      │ 0.313025    │ 0.872143    │
│ 2    │ contig_23109 │ 0.263158  │ 0.228532  │ 0.43618   │ 1.0      │ 0.508972    │ 0.987085    │
│ 3    │ contig_4493  │ 0.142857  │ 0.209751  │ 2.13589   │ 1.0      │ 0.143886    │ 0.51061     │
│ 4    │ contig_10742 │ 0.0       │ 0.0       │ 0.0       │ 0.0      │ missing     │ missing     │
│ 5    │ contig_14898 │ 0.0952381 │ 0.0907029 │ 0.0525    │ 1.0      │ 0.818769    │ 0.987085    │
│ 6    │ contig_8483  │ 0.0       │ 0.0       │ 0.0       │ 0.0      │ missing     │ missing     │
⋮
│ 2207 │ contig_18959 │ 0.52381   │ 0.427438  │ 1.06752   │ 1.0      │ 0.301507    │ 0.856472    │
│ 2208 │ contig_43517 │ 0.0952381 │ 0.0907029 │ 0.0525    │ 1.0      │ 0.818769    │ 0.987085    │
│ 2209 │ contig_27356 │ 0.047619  │ 0.0464853 │ 0.0124926 │ 1.0      │ 0.911006    │ 0.987085    │
│ 2210 │ contig_475   │ 0.285714  │ 0.30839   │ 0.113538  │ 1.0      │ 0.736152    │ 0.987085    │
│ 2211 │ contig_19384 │ 0.0952381 │ 0.176871  │ 21.0582   │ 3.0      │ 0.000102388 │ 0.000625326 │
│ 2212 │ contig_22368 │ 0.047619  │ 0.0464853 │ 0.0124926 │ 1.0      │ 0.911006    │ 0.987085    │
│ 2213 │ contig_2784  │ 0.0       │ 0.0       │ 0.0       │ 0.0      │ missing     │ missing     │, Georgia = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²         │ DF       │ P           │ Pcorr       │
│      │ String       │ Float64⍰  │ Float64⍰  │ Float64⍰   │ Float64⍰ │ Float64⍰    │ Float64⍰    │
├──────┼──────────────┼───────────┼───────────┼────────────┼──────────┼─────────────┼─────────────┤
│ 1    │ contig_35208 │ 0.4       │ 0.389467  │ 0.276208   │ 3.0      │ 0.964439    │ 1.0         │
│ 2    │ contig_23109 │ 0.28125   │ 0.304688  │ 0.378698   │ 1.0      │ 0.538301    │ 0.987085    │
│ 3    │ contig_4493  │ 0.234375  │ 0.230347  │ 0.0195733  │ 1.0      │ 0.888735    │ 0.987085    │
│ 4    │ contig_10742 │ 0.046875  │ 0.0457764 │ 0.036864   │ 1.0      │ 0.847742    │ 0.987085    │
│ 5    │ contig_14898 │ 0.292308  │ 0.249586  │ 1.90447    │ 1.0      │ 0.167579    │ 0.57153     │
│ 6    │ contig_8483  │ 0.0153846 │ 0.0152663 │ 0.00390602 │ 1.0      │ 0.950166    │ 0.996029    │
⋮
│ 2207 │ contig_18959 │ 0.538462  │ 0.461657  │ 1.79908    │ 1.0      │ 0.179824    │ 0.59782     │
│ 2208 │ contig_43517 │ 0.123077  │ 0.223432  │ 65.3201    │ 3.0      │ 4.28546e-14 │ 1.01329e-12 │
│ 2209 │ contig_27356 │ 0.0307692 │ 0.0601183 │ 65.0164    │ 3.0      │ 4.9738e-14  │ 1.01329e-12 │
│ 2210 │ contig_475   │ 0.369231  │ 0.400473  │ 0.395604   │ 1.0      │ 0.529368    │ 0.987085    │
│ 2211 │ contig_19384 │ 0.0769231 │ 0.130533  │ 68.8823    │ 3.0      │ 7.43849e-15 │ 1.01329e-12 │
│ 2212 │ contig_22368 │ 0.0923077 │ 0.116923  │ 65.1572    │ 3.0      │ 4.64073e-14 │ 1.01329e-12 │
│ 2213 │ contig_2784  │ 0.0461538 │ 0.103314  │ 65.0387    │ 3.0      │ 4.91829e-14 │ 1.01329e-12 │, South Carolina = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²        │ DF       │ P           │ Pcorr      │
│      │ String       │ Float64⍰  │ Float64⍰  │ Float64⍰  │ Float64⍰ │ Float64⍰    │ Float64⍰   │
├──────┼──────────────┼───────────┼───────────┼───────────┼──────────┼─────────────┼────────────┤
│ 1    │ contig_35208 │ 0.4       │ 0.375     │ 0.0888889 │ 1.0      │ 0.765594    │ 0.987085   │
│ 2    │ contig_23109 │ 0.117647  │ 0.110727  │ 0.0664063 │ 1.0      │ 0.796643    │ 0.987085   │
│ 3    │ contig_4493  │ 0.277778  │ 0.313272  │ 0.231066  │ 1.0      │ 0.630735    │ 0.987085   │
│ 4    │ contig_10742 │ 0.0526316 │ 0.0512465 │ 0.0138787 │ 1.0      │ 0.90622     │ 0.987085   │
│ 5    │ contig_14898 │ 0.1       │ 0.095     │ 0.0554017 │ 1.0      │ 0.813917    │ 0.987085   │
│ 6    │ contig_8483  │ 0.0       │ 0.0       │ 0.0       │ 0.0      │ missing     │ missing    │
⋮
│ 2207 │ contig_18959 │ 0.35      │ 0.43875   │ 0.818338  │ 1.0      │ 0.365667    │ 0.948296   │
│ 2208 │ contig_43517 │ 0.0       │ 0.18      │ 20.0      │ 1.0      │ 7.74422e-6  │ 5.90346e-5 │
│ 2209 │ contig_27356 │ 0.1       │ 0.185     │ 20.0617   │ 3.0      │ 0.000164815 │ 0.00078662 │
│ 2210 │ contig_475   │ 0.3       │ 0.255     │ 0.622837  │ 1.0      │ 0.429995    │ 0.987085   │
│ 2211 │ contig_19384 │ 0.0       │ 0.255     │ 20.0      │ 1.0      │ 7.74422e-6  │ 5.90346e-5 │
│ 2212 │ contig_22368 │ 0.05      │ 0.22375   │ 20.0163   │ 3.0      │ 0.000168425 │ 0.00078662 │
│ 2213 │ contig_2784  │ 0.05      │ 0.14125   │ 20.0146   │ 3.0      │ 0.000168563 │ 0.00078662 │, Florida Keys = 2213×7 DataFrame
│ Row  │ locus        │ het_obs  │ het_exp  │ χ²        │ DF       │ P        │ Pcorr    │
│      │ String       │ Float64⍰ │ Float64⍰ │ Float64⍰  │ Float64⍰ │ Float64⍰ │ Float64⍰ │
├──────┼──────────────┼──────────┼──────────┼───────────┼──────────┼──────────┼──────────┤
│ 1    │ contig_35208 │ 0.45     │ 0.41125  │ 0.519382  │ 3.0      │ 0.914613 │ 0.987085 │
│ 2    │ contig_23109 │ 0.4      │ 0.42     │ 0.0453515 │ 1.0      │ 0.831359 │ 0.987085 │
│ 3    │ contig_4493  │ 0.25     │ 0.28875  │ 0.360188  │ 1.0      │ 0.548402 │ 0.987085 │
│ 4    │ contig_10742 │ 0.2      │ 0.18     │ 0.246914  │ 1.0      │ 0.619257 │ 0.987085 │
│ 5    │ contig_14898 │ 0.25     │ 0.21875  │ 0.408163  │ 1.0      │ 0.522903 │ 0.987085 │
│ 6    │ contig_8483  │ 0.05     │ 0.04875  │ 0.0131492 │ 1.0      │ 0.908707 │ 0.987085 │
⋮
│ 2207 │ contig_18959 │ 0.45     │ 0.39875  │ 0.330382  │ 1.0      │ 0.565434 │ 0.987085 │
│ 2208 │ contig_43517 │ 0.05     │ 0.04875  │ 0.0131492 │ 1.0      │ 0.908707 │ 0.987085 │
│ 2209 │ contig_27356 │ 0.0      │ 0.0      │ 0.0       │ 0.0      │ missing  │ missing  │
│ 2210 │ contig_475   │ 0.4      │ 0.375    │ 0.0888889 │ 1.0      │ 0.765594 │ 0.987085 │
│ 2211 │ contig_19384 │ 0.15     │ 0.13875  │ 0.131483  │ 1.0      │ 0.7169   │ 0.987085 │
│ 2212 │ contig_22368 │ 0.15     │ 0.13875  │ 0.131483  │ 1.0      │ 0.7169   │ 0.987085 │
│ 2213 │ contig_2784  │ 0.0      │ 0.0      │ 0.0       │ 0.0      │ missing  │ missing  │, Mideast Gulf = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp  │ χ²        │ DF       │ P           │ Pcorr      │
│      │ String       │ Float64⍰  │ Float64⍰ │ Float64⍰  │ Float64⍰ │ Float64⍰    │ Float64⍰   │
├──────┼──────────────┼───────────┼──────────┼───────────┼──────────┼─────────────┼────────────┤
│ 1    │ contig_35208 │ 0.5       │ 0.436224 │ 0.598475  │ 1.0      │ 0.439161    │ 0.987085   │
│ 2    │ contig_23109 │ 0.214286  │ 0.191327 │ 0.4032    │ 1.0      │ 0.525441    │ 0.987085   │
│ 3    │ contig_4493  │ 0.178571  │ 0.21875  │ 0.944606  │ 1.0      │ 0.331096    │ 0.898252   │
│ 4    │ contig_10742 │ 0.107143  │ 0.101403 │ 0.0897116 │ 1.0      │ 0.764544    │ 0.987085   │
│ 5    │ contig_14898 │ 0.357143  │ 0.336735 │ 0.102847  │ 1.0      │ 0.74844     │ 0.987085   │
│ 6    │ contig_8483  │ 0.0       │ 0.0      │ 0.0       │ 0.0      │ missing     │ missing    │
⋮
│ 2207 │ contig_18959 │ 0.535714  │ 0.484056 │ 0.318893  │ 1.0      │ 0.572274    │ 0.987085   │
│ 2208 │ contig_43517 │ 0.178571  │ 0.162628 │ 0.269127  │ 1.0      │ 0.603918    │ 0.987085   │
│ 2209 │ contig_27356 │ 0.0357143 │ 0.101403 │ 11.75     │ 1.0      │ 0.000608429 │ 0.00278639 │
│ 2210 │ contig_475   │ 0.285714  │ 0.336735 │ 0.642792  │ 1.0      │ 0.422702    │ 0.987085   │
│ 2211 │ contig_19384 │ 0.0       │ 0.0      │ 0.0       │ 0.0      │ missing     │ missing    │
│ 2212 │ contig_22368 │ 0.107143  │ 0.101403 │ 0.0897116 │ 1.0      │ 0.764544    │ 0.987085   │
│ 2213 │ contig_2784  │ 0.0714286 │ 0.135204 │ 28.0414   │ 3.0      │ 3.56005e-6  │ 3.64412e-5 │, Northeast Gulf = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²        │ DF       │ P          │ Pcorr      │
│      │ String       │ Float64⍰  │ Float64⍰  │ Float64⍰  │ Float64⍰ │ Float64⍰   │ Float64⍰   │
├──────┼──────────────┼───────────┼───────────┼───────────┼──────────┼────────────┼────────────┤
│ 1    │ contig_35208 │ 0.428571  │ 0.375     │ 0.571429  │ 1.0      │ 0.449692   │ 0.987085   │
│ 2    │ contig_23109 │ 0.321429  │ 0.392219  │ 0.912124  │ 1.0      │ 0.339552   │ 0.903179   │
│ 3    │ contig_4493  │ 0.214286  │ 0.191327  │ 0.4032    │ 1.0      │ 0.525441   │ 0.987085   │
│ 4    │ contig_10742 │ 0.0357143 │ 0.0350765 │ 0.0092562 │ 1.0      │ 0.923354   │ 0.987085   │
│ 5    │ contig_14898 │ 0.25      │ 0.21875   │ 0.571429  │ 1.0      │ 0.449692   │ 0.987085   │
│ 6    │ contig_8483  │ 0.0357143 │ 0.0350765 │ 0.0092562 │ 1.0      │ 0.923354   │ 0.987085   │
⋮
│ 2207 │ contig_18959 │ 0.285714  │ 0.336735  │ 0.642792  │ 1.0      │ 0.422702   │ 0.987085   │
│ 2208 │ contig_43517 │ 0.107143  │ 0.101403  │ 0.0897116 │ 1.0      │ 0.764544   │ 0.987085   │
│ 2209 │ contig_27356 │ 0.142857  │ 0.132653  │ 0.16568   │ 1.0      │ 0.68398    │ 0.987085   │
│ 2210 │ contig_475   │ 0.464286  │ 0.448342  │ 0.0354101 │ 1.0      │ 0.850739   │ 0.987085   │
│ 2211 │ contig_19384 │ 0.0357143 │ 0.0350765 │ 0.0092562 │ 1.0      │ 0.923354   │ 0.987085   │
│ 2212 │ contig_22368 │ 0.0357143 │ 0.0350765 │ 0.0092562 │ 1.0      │ 0.923354   │ 0.987085   │
│ 2213 │ contig_2784  │ 0.0       │ 0.0688776 │ 28.0      │ 1.0      │ 1.21315e-7 │ 1.87007e-6 │, Southeast Gulf = 2213×7 DataFrame
│ Row  │ locus        │ het_obs   │ het_exp   │ χ²         │ DF       │ P          │ Pcorr      │
│      │ String       │ Float64⍰  │ Float64⍰  │ Float64⍰   │ Float64⍰ │ Float64⍰   │ Float64⍰   │
├──────┼──────────────┼───────────┼───────────┼────────────┼──────────┼────────────┼────────────┤
│ 1    │ contig_35208 │ 0.433333  │ 0.375     │ 0.725926   │ 1.0      │ 0.394207   │ 0.980648   │
│ 2    │ contig_23109 │ 0.2       │ 0.231111  │ 0.543639   │ 1.0      │ 0.460928   │ 0.987085   │
│ 3    │ contig_4493  │ 0.133333  │ 0.18      │ 2.01646    │ 1.0      │ 0.155601   │ 0.540012   │
│ 4    │ contig_10742 │ 0.0666667 │ 0.0644444 │ 0.0356718  │ 1.0      │ 0.850195   │ 0.987085   │
│ 5    │ contig_14898 │ 0.2       │ 0.18      │ 0.37037    │ 1.0      │ 0.542802   │ 0.987085   │
│ 6    │ contig_8483  │ 0.0333333 │ 0.0327778 │ 0.00861821 │ 1.0      │ 0.926035   │ 0.987085   │
⋮
│ 2207 │ contig_18959 │ 0.466667  │ 0.42      │ 0.37037    │ 1.0      │ 0.542802   │ 0.987085   │
│ 2208 │ contig_43517 │ 0.1       │ 0.095     │ 0.0831025  │ 1.0      │ 0.773136   │ 0.987085   │
│ 2209 │ contig_27356 │ 0.133333  │ 0.124444  │ 0.153061   │ 1.0      │ 0.695627   │ 0.987085   │
│ 2210 │ contig_475   │ 0.433333  │ 0.375     │ 0.725926   │ 1.0      │ 0.394207   │ 0.980648   │
│ 2211 │ contig_19384 │ 0.0666667 │ 0.0644444 │ 0.0356718  │ 1.0      │ 0.850195   │ 0.987085   │
│ 2212 │ contig_22368 │ 0.133333  │ 0.124444  │ 0.153061   │ 1.0      │ 0.695627   │ 0.987085   │
│ 2213 │ contig_2784  │ 0.0       │ 0.124444  │ 30.0       │ 1.0      │ 4.32046e-8 │ 7.21496e-7 │)
```



### pro tip 

If using `by_pop = true`, there may be a very long output which you may want to suppress by ending the command with a semicolon `;`.   The function returns a NamedTuple of DataFrames, meaning you can index it with a dot `.` operator. If you assign the command's output to a variable, such as `hardy = hwe_test(gulfsharks(), by_pop = true) ;`, then you can index the output with the population names, such as `hardy.CapeCanaveral` or `hardy.Georgia`.  And since those objects are DataFrames, you can continue using the dot `.` operator to index them further by column, such as `hardy.Georgia.Pcorr` to see only the corrected P values. Convenience!



!!!info "indexing the Chiq-sq column"
    If you wish to use the dot operator to index the chi-squared values, you will need to use the unicode characters in Julia to do so, b/c the column is literally named `χ²`. To generate those characters, type in `\Chi` + press `TAB` + type in `\^2` + press `TAB` without spaces and it will magically appear. Written out in more explicit English (we really want you to get it!), it's a backslash `\`, the word `Chi` with a capital C, the TAB key on your keyboard (and you'll notice it's already changed it to the letter χ), another backslash `\`, a caret `^`, the number `2`, then the TAB key again. 
    
    All in all, you'll be doing this: 
    
    `\ChiTAB\^2TAB`. It's actually a lot easier than it looks.
