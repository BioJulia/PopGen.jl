---
id: generalutils
title: GeneralUtils.jl
sidebar_label: GeneralUtils.jl
---
## PopGenCore.jl/src/Utils/GeneralUtils.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪🔵 Base.copy
```julia
Base.copy(data::PopData)
```

----
### 🟪🔵 Base.size
```julia
Base.size(data::PopData)
```

----
### 🟪🔵 Base.sort
```julia
Base.sort(x::NTuple{N,T}) where N where T <: Signed 
```
Sort the integers within a Tuple and return the sorted Tuple.

----
### 🟪 convertcoord
```julia
convertcoord(coordinate::String)
```
Takes non-decimal-degree format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of coordinate strings to convert them.
##### Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)

**Example**
```
julia> convertcoord("-41 31.52")
-41.5253f0
julia> convertcoord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```

----
### 🟪🔵 dropmonomorphic
```julia
dropmonomorphic(data::PopData; silent::Bool = false)
```
Return a `PopData` object omitting any monomorphic loci. Will inform you which
loci were removed.

----
### 🟪🔵 dropmonomorphic!
```julia
dropmonomorphic!(data::PopData; silent::Bool = false)
```
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which
loci were removed.

----
### 🟪🔵 dropmultiallelic(data::PopData)
```julia
dropmultiallelic
```
Return a `PopData` object omitting loci that are not biallelic.

----
### 🟪🔵 dropmultiallelic!(data::PopData)
```julia
dropmultiallelic!
```
Edit a `PopData` object in place, removing loci that are not biallelic.

----
### ❗truncatepath
```julia
truncatepath(text::String)
```