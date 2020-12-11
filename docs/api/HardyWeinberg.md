---
id: hardyweinberg
title: HardyWeinberg.jl
sidebar_label: HardyWeinberg.jl
---

### `locus_chi_sq`
```julia
locus_chi_sq(locus::T) where T <: GenotypeArray
```
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value.

----

### `hwe_test`
```julia
    hwe_test(data::PopData; by_pop::Bool = false; correction = "none")
```
Calculate chi-squared test of HWE for each locus and returns observed and
expected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use `by_pop = true` to perform this separately for each population (default: by_pop = false) and return a NamedTuple of DataFrames with the names corresponding to the population names. Use `correction =` to specify a P-value adjustment method for multiple testing.

**correction methods (case insensitive)**
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

**Example**
```julia
hwe_test(@gulfsharks, correction = "bh")
hwe_test(@gulfsharks, by_pop = true, correction = "bh")
```
