Like most Julia packages, there are a lot of functions in PopGen.jl, but only a number of them are `exported`, which are the ones you the user can casually use after calling `using PopGen`. You'll often see other Julia packages refer to these as API's, or "Application Programming Interface". If you want to contribute to PopGen.jl, then it will be super useful to know what API's are already available and save you the trouble of reinventing the wheel.



## Allele Frequencies

These are found in `AlleleFreq.jl`

### allele_freq_mini

```julia
allele_freq_mini(x::Array{Union{Missing, Tuple},1})
```

Calculate allele counts for a single locus of a `PopObj`. Returns a `Dict` of allele's and their frequencies.

```julia
allele_freq_mini(x::Array{Union{Missing, Tuple},1})
```

Calculate allele counts for a single locus of a `PopObj` split by population using `group()`. Returns a `Dict` of allele's and their frequencies.



## Genotype Frequencies

These are found in `AlleleFreq.jl`

### geno_freq

```julia
geno_freq(x::Array{Union{Missing, Tuple},1})
```

Calculate genotype frequencies of all loci in a `PopObj`. Returns a `Dict` of genotypes and their frequencies. 

```julia
geno_freq(x::SubArray{Union{Missing, Tuple},1})
```

Calculate genotype frequencies of all loci in `PopObj` split by population using `group()`. Returns a `Dict` of genotypes and their frequencies.



## Observed Heterozygosity

These are found in `HardyWeinberg.jl`

### het_observed

```julia
het_observed(x::PopObj)
```

Calculate the observed heterozygosity for each locus in a `PopObj`. Returns an array of heterozygosity values.

### het_population_obs

```julia
het_population_obs(x::PopObj)
```

Return a `Dict` of the observed heterozygosity per population for each locus in a `PopObj`

### het_sample

```julia
het_sample(x::PopObj)
```

Calculate the observed heterozygosity for each individual in a `PopObj`. Returns an array of heterozygosity values.



## Expected Heterozygosity

These are found in `HardyWeinberg.jl`

### het_expected

```julia
het_expected(x::PopObj)
```

Calculate the expected heterozygosity for each locus in a `PopObj`. Returns an array of heterozygosity values.

### het_population_exp

```julia
het_population_exp(x::PopObj)
```

Return a `Dict` of the expected heterozygosity per population for each locus in a `PopObj`.



## Chi Squared

Found in `HardyWeinberg.jl`

### locus_chi_sq

```julia
locus_chi_sq(locus::Array{Union{Missing, Tuple},1})
```

Calculate the chi square statistic and p-value for a locus. Returns a tuple with chi-square statistic, degrees of freedom, and p-value.

## Multiple Testing

Found in `HardyWeinberg.jl`

### multitest_missing

```julia
multitest_missing(pvals::Array{Float64,1}, correction::String)
```

Modification to `MultipleTesting.adjust` to include `missing` values. Missing values are first removed from the array, the appropriate correction made, then missing values are re-added to the array at their original positions. See MultipleTesting.jl docs for full more detailed information.

Example:

`multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")`

`correction` methods (case insensitive):

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