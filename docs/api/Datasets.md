---
id: datasets
title: Datasets.jl
sidebar_label: Datasets.jl
---

### `dataset`
```julia
dataset(::String)
```
Load an example dataset from either `"gulfsharks"` (SNP) or `"nancycats"` (microsatellite). Can also use `"sharks"` and `"cats"`
as shorthands. Use `?nancycats` and `?gulfsharks` to learn more about
these datasets.

**Example**
```
ncats = dataset("nancycats")
gsharks = dataset("sharks")
```

----

### `nancycats`
```julia
@nancycats
```
Returns a `PopObj` of corresponding "nancycats" dataset as featured in
the R package `adegenet`. This is microsatellite data of 9 loci in 237
individuals across 17 populations.

**Example**
```
ncats = @nancycats
```

----

### `gulfsharks`
```julia
@gulfsharks
```
Returns a `PopObj` corresponding the Blacknose shark dataset as used in
Dimens et al. 2019. This is a mid-sized SNP dataset of 2213 SNPs in 212
individuals, across 7 populations.

**Example**
```julia
sharks = @gulfsharks
```