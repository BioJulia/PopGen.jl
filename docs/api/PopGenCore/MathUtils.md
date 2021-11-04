---
id: mathutils
title: MathUtils.jl
sidebar_label: MathUtils.jl
---
## PopGenCore.jl/src/Utils/MathUtils.jl
📦  => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪 countnonzeros
```julia
countnonzeros(x::AbstractVector{T}) where T<:Real
```
Return the number of non-zero values in a vector

----
### 🟪 reciprocal
```julia
reciprocal(num::T) where T <: Signed
```
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.

----
### 🟪 reciprocalsum
```julia
reciprocalsum(x::AbstractVector{T}) where T<:Real
```
Return the sum of the reciprocal values of `x`, skipping the `Inf` values
resulting from divide-by-zero errors.
