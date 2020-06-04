---
id: utils
title: Utils.jl
sidebar_label: Utils.jl
---

### `convert_coord`
```julia
convert_coord(coordinate::string)
```
Takes non-decimal-degree format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of coordinate strings to convert them.
#### Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
**Example**
```
julia> convert_coord("-41 31.52")
-41.5253f0
julia> convert_coord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```

### `nonmissing`
```julia
nonmissing(vec::T) where T<:AbstractArray
```
Convenience function to count the number of non-`missing` values in a vector.


### `reciprocal`
```julia
reciprocal(num::T) where T <: Signed
```
Returns the reciprocal (1/number) of a number. Will return `0` when 
the number is `0` instead of returning `Inf`.

### `multitest_missing`
```julia
multitest_missing(pvals::Array{Float64,1}, correction::String)
```
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. Missing values are first removed from the array, the appropriate
correction made, then missing values are re-added to the array at their original
positions. See MultipleTesting.jl docs for full more detailed information.

**Example**
```julia
multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")`
```
#### `correction` methods (case insensitive)
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"` : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
