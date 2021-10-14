---
id: conditionals
title: Conditionals and Logic
sidebar_label: Conditionals
---

Included in PopGen.jl are some functions to help discriminate your data a bit more. Like all conditionals, these return `true` or `false` depending on the test. 

## Homozygosity
```julia
ishom(locus::Genotype)
ishom(locus::GenoArray)
```
This will return `true` if a genotype is homozygous. The `GenoaArray` version
just broadcasts it across all the genotypes in an array, returning a vector
of `true` or `false`. Returns `missing` if the genotype is `missing`.

**Example**
```julia
julia> cats = @nancycats ;
julia> subset = cats.genodata.genotype[1:10]
10-element Vector{Union{Missing, Tuple{Int16, Int16}}}:
 missing
 missing
 (135, 143)
 (133, 135)
 (133, 135)
 (135, 143)
 (135, 135)
 (135, 143)
 (137, 143)
 (135, 135)

julia> ishom(subset[3])
false

julia> ishom(subset)
10-element Vector{Union{Missing, Bool}}:
      missing
      missing
 false
 false
 false
 false
  true
 false
 false
  true

```
:::note using skipmissing
If you want to avoid `missing` genotypes, you can use `skipmissing` to ignore them. This also works for `ishet`.
```julia
julia> ishom(skipmissing(subset))
8-element Vector{Bool}:
 0
 0
 0
 0
 1
 0
 0
 1
```
:::

Another option is to check if a genotype is homozygous for a specific allele. To
do that, we exploit Julia's multiple dispatch and use `ishom` again, but with
different arguments.
```julia
ishom(geno::Genotype, allele::Signed)
ishom(genos::GenoArray, allele::Signed)
```
This will return `true` if the `geno` (or `genos`) is/are homozygous for the specified `allele`. Notice that when we query a genotype that doesn't contain that allele, it returns `false`.

**Example**
```julia
julia> ishom(subset[3], 135)
false

julia> ishom(subset[10], 135)
true

julia> ishom(subset[9], 135)
false

julia> ishom(subset, 135)
10-element Vector{Union{Missing, Bool}}:
      missing
      missing
 false
 false
 false
 false
  true
 false
 false
  true
```


## Heterozygosity
```julia
ishet(locus::Genotype)
ishet(locus::GenoArray)
```
This is the exact opposite of `ishom`, returning `true` if the genotype (or genotypes) is/are heterozygous. Returns `missing` if the genotype is `missing`.

**Example**
```julia
julia> cats = @nancycats ;
julia> subset = cats.genodata.genotype[1:10]
10-element Vector{Union{Missing, Tuple{Int16, Int16}}}:
 missing
 missing
 (135, 143)
 (133, 135)
 (133, 135)
 (135, 143)
 (135, 135)
 (135, 143)
 (137, 143)
 (135, 135)

julia> ishet(subset[3])
true

julia> ishet(subset)
10-element Vector{Union{Missing, Bool}}:
      missing
      missing
  true
  true
  true
  true
 false
  true
  true
 false
 ```

We likewise have the option to check if a locus is heterozygous for a specific
allele. To do that, we again exploit Julia's multiple dispatch and use `ishet`, 
but with different arguments.

```julia
ishet(geno::Genotype, allele::Signed)
ishet(genos::GenoArray, allele::Signed)
```
This will return `true` if the `geno` (or `genos`) is/are heterozygous for the specified `allele`. Notice that when we query a genotype that doesn't contain that allele, it returns `false`.

**Example**
```julia
julia> ishet(subset[3], 135)
true

julia> ishet(subset[10], 135)
false

julia> ishet(subset[9], 135)
true

julia> ishom(subset, 135)
10-element Vector{Union{Missing, Bool}}:
      missing
      missing
  true
  true
  true
  true
 false
  true
 false
 false
```

## Biallelic data
Some analyses are restricted to work exclusively on biallelic data (e.g. Hudson pairwise FST), so it may help to know if things are biallelic.

```julia
isbiallelic(data::GenoArray)
```
Returns `true` if the `GenoArray` is biallelic, `false` if not.

```julia
isbiallelic(data::PopData)
```
Returns `true` all the loci in the `PopData` are biallelic, `false` if not.

**Example**
```julia
julia> sharks = @gulfsharks ;

julia> isbiallelic(sharks)
false

julia> drop_multiallelic!(sharks)
[ Info: Removing 258 multialleic loci

julia> isbiallelic(sharks)
true
```

