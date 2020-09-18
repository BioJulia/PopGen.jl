---
id: popgensims_cross
title: Cross.jl
sidebar_label: Cross.jl
---

### `sample_genotype`
```julia
sample_genotype(geno::T, n_alleles::Int) where T<:Genotype
```
```julia
sample_genotype(geno::Missing, n_alleles::Int)
```
----

### `haploid_cross!`
```julia
haploid_cross!(data::DataFrame, p1::T, p2::T; n::Int) where T <: GenoArray
```

----

### `polyploid_cross!`
```julia
polyploid_cross!(data::DataFrame, p1::T, p2::T; n::Int, ploidy::Int) where T <: GenoArray
 ```

----
    
### `cross`
```julia
cross(data::PopData, parent1::String, parent2::String; n::Int = 100, generation::String = "F1")
```
Simulate a breeding cross between individuals `parent1` and `parent2` from a `PopData` object.
Returns PopData consisting of `n` offspring resulting from the cross.
#### Keyword Arguments
- `n` : Integer of number of offspring to generate (default: `100`)
- `generation` : A string to assign `population` identity to the offspring (default: `"F1"`)



```julia
cross(parent_1::Pair, parent_2::Pair, n::Int = 100, generation::String = "F1")
```
Simulate a breeding cross between individuals `parent` and `parent2` from two different `PopData` objects.
Returns PopData consisting of `n` offspring resulting from the cross. `parent_1_data` and `parent_2_data` 
are positional arguments, therefore they must be written without keywords and in the order of parents 1, parent 2. 
#### Keyword Arguments
- `parent_1` : Pair of `PopData => "Parent1Name"`
- `parent_2` : Pair of `PopData => "Parent1Name"`
- `n` : Integer of number of offspring to generate (default: `100`)
- `generation` : A string to assign `population` identity to the offspring (default: `"F1"`)
