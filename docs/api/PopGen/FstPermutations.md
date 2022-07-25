---
id: fstpermutations
title: FstPermutations.jl
sidebar_label: FstPermutations.jl
---
## PopGen.jl/src/FStatistics/FstPermutations.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 📦 _fst_permutate
```julia
_fst_permute(population_1::T, population_2::T) where T<:AbstractMatrix
```
Returns two matrices with rows (samples) shuffled between them. Respects the
number of rows of the original matrices (i.e. population sizes).

----

### 📦 _fst_permutation
```julia
_fst_permution(data::PopData, method::Function, iterations::Int64)
```
Returns a `PairwiseFST` object containing a dataframe of Pairwise FST calculations. The contained
dataframe has FST values below the diagonal and P values above it. This method is used internally
and wrapped by the public API provided in `pairwisefst()`.
