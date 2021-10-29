---
id: conditionals
title: Conditionals.jl
sidebar_label: Conditionals.jl
---
## PopGenCore.jl/src/Conditionals.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪🔵 isbiallelic
```julia
isbiallelic(data::GenoArray)
```
Returns `true` if the `GenoArray` is biallelic, `false` if not.

----

```julia
isbiallelic(data::DataFrame)
```
----

```julia
isbiallelic(data::PopData)
```
Returns `true` all the loci in the `PopData` are biallelic, `false` if not.

----
### 🟪🔵 ishom
```
ishom(locus::T) where T <: GenoArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't (or missing). For calculations, we recommend using `_ishom()`,
which returns `missing` if the genotype is `missing`. The vector version
simply maps the function over the elements.

----

```julia
ishom(locus::Genotype, allele::Signed)
ishom(loci::GenoArray, allele::Signed)
ishom(geno::Missing, allele::Signed)
```
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele`.

### 🟪 _ishom
```julia
_ishom(locus::T) where T <: GenoArray
_ishom(locus::Genotype)
_ishom(locus::Missing)
```
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele` and
`missing` if the genotype is `missing`.

### 🟪🔵 ishet
```julia
ishet(locus::T) where T <: GenoArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if
it is, `false` if it isn't. The vector version simply broadcasts the function over the
elements. Under the hood, this function is simply `!ishom`.
function ishet(locus::Genotype)


```julia
ishet(locus::Genotype, allele::Signed)
ishet(loci::GenoArray, allele::Signed)
```
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele`. 


### 🟪 _ishet
```
_ishet(locus::T) where T <: GenoArray
_ishet(locus::Genotype)
_ishet(locus::Missing)
```
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele` and
`missing` if the genotype is `missing`.