---
id: permutations
title: Permutations.jl
sidebar_label: Permutations.jl
---

### `permute_loci!`
```julia
permute_loci!(data::PopData)
```
Edits `PopData` in place with loci permuted across populations within
the `.loci` dataframe.

----

### `permute_samples!`
```julia
permute_samples!(data::PopData; meta::Bool = false)
```
Edits `PopData` in place with samples permuted across populations within
the `.loci` dataframe. Since performance is important for many permutations,
the default is to only edit the `.loci` table in place; use `meta = true`
if you also require the `.meta` dataframe edited in place.

----

### `permute_genotypes!`
```julia
permute_genotypes!(data::PopData; by::String = "locus")
```
Edits `PopData` in place with genotypes permuted across individuals within
the `.loci` dataframe. Use `by = "population"` (or `"pop"`) to permute genotypes
within populations.

----

### `permute_alleles!`
```julia
permute_alleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
```
Edits `PopData` in place with alleles permuted and reconstructed into genotypes
for each locus within the `.loci` dataframe. Use `by = "population"` (or `"pop"`)
to permute alleles within populations. If `ploidy` is not provided (default `ploidy = nothing`),
then ploidy will be identified from the PopData. If performance is important,
it would be best to identify ploidy in advance and set it to a specific integer.
