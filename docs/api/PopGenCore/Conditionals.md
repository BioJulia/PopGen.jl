---
id: conditionals
title: Conditionals.jl
sidebar_label: Conditionals.jl
---

### `ishom`
```julia
ishom(locus::T) where T <: GenotypeArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if it is, `false` if it isn't, and `missing` if it's `missing`. The vector version simply broadcasts the function over the elements.

```julia
ishom(locus::Genotype, allele::Signed)
ishom(loci::GenoArray, allele::Signed)
```
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele`.

----

### `ishet`
```julia
ishet(locus::T) where T <: GenotypeArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if it is, `false` if it isn't. The vector version simply broadcasts the function over the elements. Under the hood, this function is simply `!ishom`.

```julia
ishet(locus::Genotype, allele::Signed)
ishet(loci::GenoArray, allele::Signed)
```
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele`. 

----

### `isbiallelic`
```julia
isbiallelic(data::T) where T<:GenoArray
```
Returns `true` if the `GenoArray` is biallelic, `false` if not.

```julia
isbiallelic(data::PopData)
```
Returns `true` all the loci in the `PopData` are biallelic, `false` if not.

----