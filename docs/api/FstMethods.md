---
id: fstmethods
title: FstMethods.jl
sidebar_label: FstMethods.jl
---

### `_pairwise_Hudson`
```
_pairwise_Hudson(data::PopData)
```

----
### `hudson_fst`
```julia
hudson_fst(population_1::T, population_2::T) where T<:AbstractMatrix
```
----

### `_hudson_fst`
```julia
_hudson_fst(pop1::T, pop2::T) where T<:GenoArray
```
This is a helper function to perform the math on a single locus

-----

### `nei_fst`
```julia
nei_fst(population_1::T, population_2::T) where T<:AbstractMatrix
```

----

### `_pairwise_Nei`
```julia
_pairwise_Nei(data::PopData)
```

----

### `weircockerham_fst`
```julia
weircockerham_fst(population_1::T, population_2::T) where T<:AbstractMatrix
```
----

### `_pairwise_WeirCockerham`
```julia
_pairwise_WeirCockerham(data::PopData)
```