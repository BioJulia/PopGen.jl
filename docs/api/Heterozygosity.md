---
id: heterozygosity
title: Heterozygosity.jl
sidebar_label: Heterozygosity.jl
---

### `ishom`
```julia
ishom(locus::T) where T <: GenotypeArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if it is, `false` if it isn't, and `missing` if it's `missing`. The vector version simply broadcasts the function over the elements.

### `ishet`
```julia
ishet(locus::T) where T <: GenotypeArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if it is, `false` if it isn't. The vector version simply broadcasts the function over the elements. Under the hood, this function is simply `!ishom`.

### `hetero_o`
```julia
hetero_o(data::T) where T <: GenotypeArray
```
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined as genotypes returning `true` for `ishet()`. This is numerically feasible because `true` values are mathematically represented as `1`, whereas `false` are represented as `0`.

### `hetero_e`
```julia
hetero_e(allele_freqs::Vector{T}) where T <: GenotypeArray
```
Returns the expected heterozygosity of an array of genotypes, calculated as 1 - sum of the squared allele frequencies.

### `heterozygosity`
```julia
heterozygosity(data::PopData, by::String = "locus")
```
Calculate observed and expected heterozygosity in a `PopData` object. For loci, heterozygosity is calculated in the Nei fashion, such that heterozygosity is calculated as the average over heterozygosity per locus per population.

**Modes**
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population

**Example**
```julia
heterozygosity(nancycats(), "population" )
```

### `het_sample`
```julia
    het_sample(data::PopData, individual::String)
```
Calculate the observed heterozygosity for an individual in a `PopData` object. Returns an array of heterozygosity values.
