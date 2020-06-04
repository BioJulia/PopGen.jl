---
id: allelefreq
title: AlleleFreq.jl
sidebar_label: AlleleFreq.jl
---

### `alleles`
```julia
alleles(locus::T) where T<:GenotypeArray
```
Return an array of all the non-missing alleles present in a locus.

### `unique_alleles`
```julia
unique_alleles(locus::T) where T<:GenotypeArray
```
Return an array of all the unique non-missing alleles of a locus.

### `allele_freq`
```julia
allele_freq(locus::T) where T<:GenotypeArray
```
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object.

### `allele_freq`
```julia
allele_freq(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of allele frequencies
of that locus per population.

**Example**
```julia
cats = nancycats()
allele_freq(cats, "fca8")
allele_freq(cats, "fca8", population = true)
```

### `allele_freq_vec`
```julia
allele_freq_vec(locus::T) where T<:GenotypeArray
```
Return a Vector of allele frequencies of a single locus in a `PopData` object. Similar to `allele_freq()`, except it returns only the frequencies, without the allele names, meaning they can be in any order. This can be useful for getting the expected genotype frequencies.

### `geno_count_observed`
```julia
geno_count_observed(locus::T) where T<:GenotypeArray
```
Return a `Dict` of genotype counts of a single locus in a
`PopData` object.

### `geno_count_expected`
```julia
geno_count_expected(locus::T) where T<:GenotypeArray
```
Return a `Dict` of the expected genotype counts of a single locus in a
`PopData` object. Expected counts are calculated as the product of observed allele frequencies multiplied by the number of non-missing genotypes.

### `geno_freq`
```julia
geno_freq(locus::T) where T<:GenotypeArray
`PopData` object.
```
Return a `Dict` of genotype frequencies of a single locus in a

### `geno_freq`
```julia
geno_freq(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.

```julia
cats = nancycats()
geno_freq(cats, "fca8")
geno_freq(cats, "fca8", population = true)
```

### `geno_freq`
```julia
geno_freq_expected(locus::T) where T<:GenotypeArray
```
Return a `Dict` of the expected genotype frequencies of a single locus in a `PopData` object. Expected frequencies are calculated as the product of
observed allele frequencies.

### `geno_freq_expected`
```julia
geno_freq_expected(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of expected genotype frequencies of a single locus in a
`PopData` object. Use `population = true` to return a table of expected genotype frequencies of that locus per population.

**Example**
```
cats = nancycats()
geno_freq_expeced(cats, "fca8")
geno_freq_expected(cats, "fca8", population = true)
```