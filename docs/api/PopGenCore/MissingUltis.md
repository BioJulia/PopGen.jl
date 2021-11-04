---
id: missingutils
title: MissingUtils.jl
sidebar_label: MissingUtils.jl
---
## PopGenCore.jl/src/Utils/MissingUtils.jl
📦  => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl


### 🟪 nonmissing
```julia
nonmissing(vec::T) where T<:AbstractArray
```
Convenience function to count the number of non-`missing` values
in a vector.

----
```julia
nonmissing(data::PopData, locus::String)
```
Convenience function to count the number of non-`missing` samples
at a locus.

----
### 🟪 nonmissings
```julia
nonmissings(vec1::AbstractVector, vec2::AbstractVector)
```
Return a vector of indices where neither input vectors have a `missing` value, i.e. an
intersection of the indices of their non-missing elements.
