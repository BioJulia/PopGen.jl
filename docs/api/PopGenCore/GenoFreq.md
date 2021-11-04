---
id: genofreq
title: GenoFreq.jl
sidebar_label: GenoFreq.jl
---
## PopGenCore.jl/src/GenoFreq.jl
📦  => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪 genocount_observed
```julia
genocount_observed(locus::T) where T<:GenotypeArray
```
Return a `Dict` of genotype counts of a single locus in a
`PopData` object.

----

### 🟪 genocount_expected
```julia
genocount_expected(locus::T) where T<:GenotypeArray
```
Return a `Dict` of the expected genotype counts of a single locus in a
`PopData` object. Expected counts are calculated as the product of observed allele frequencies multiplied by the number of non-missing genotypes.

----

### 🟪 genofreq
```julia
genofreq(locus::T) where T<:GenotypeArray
`PopData` object.
```
Return a `Dict` of genotype frequencies of a single locus in a

----

```julia
genofreq(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.

```julia
cats = @nancycats
genofreq(cats, "fca8")
genofreq(cats, "fca8", population = true)
```

----
### 🟪 genofreq_expected
```julia
genofreq_expected(locus::T) where T<:GenotypeArray
```
Return a `Dict` of the expected genotype frequencies of a single locus in a `PopData` object. Expected frequencies are calculated as the product of
observed allele frequencies.

----

```julia
genofreq_expected(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of expected genotype frequencies of a single locus in a
`PopData` object. Use `population = true` to return a table of expected genotype frequencies of that locus per population.

**Example**
```
cats = @nancycats
genofreq_expeced(cats, "fca8")
genofreq_expected(cats, "fca8", population = true)
```