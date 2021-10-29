---
id: permutations
title: Permutations.jl
sidebar_label: Permutations.jl
---
## PopGenCore.jl/src/Permutations.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪 permuteloci!
```julia
permuteloci!(data::PopData)
```
Edits `PopData` in place with loci permuted across populations within
the `.genodata` dataframe.

----

### 🟪 permutesamples!
```julia
permutesamples!(data::PopData; meta::Bool = false)
```
Edits `PopData` in place with samples permuted across populations within
the `.genodata` dataframe. Since performance is important for many permutations,
the default is to only edit the `.genodata` table in place; use `meta = true`
if you also require the `.sampleinfo` dataframe edited in place.

----

### 🟪 permutegenotypes!
```julia
permutegenotypes!(data::PopData; by::String = "locus")
```
Edits `PopData` in place with genotypes permuted across individuals within
the `.genodata` dataframe. Use `by = "population"` (or `"pop"`) to permute genotypes
within populations.

----

### 🟪 permutealleles!
```julia
permutealleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
```
Edits `PopData` in place with alleles permuted and reconstructed into genotypes
for each locus within the `.genodata` dataframe. Use `by = "population"` (or `"pop"`)
to permute alleles within populations. If `ploidy` is not provided (default `ploidy = nothing`),
then ploidy will be identified from the PopData. If performance is important,
it would be best to identify ploidy in advance and set it to a specific integer.


### 🟪 strictshuffle
```julia
strictshuffle(x::T) where T <: AbstractArray
```
Shuffle only the non-missing values of a Vector and return a copy of the vector,
keeping the `missing` values at their original locations.
Use `strictshuffle!` to edit in-place instead of returning a copy.

----

### 🟪 strictshuffle!
```julia
strictshuffle!(x::T) where T <: AbstractArray
```
Shuffle only the non-missing values of a Vector, keeping the
`missing` values at their original locations. Use `strictshuffle`
to return a copy instead of editing in-place.
