---
id: utils
title: Utils.jl
sidebar_label: Utils.jl
---

### `alleles`
```julia
alleles(locus::T; miss::Bool = false) where T<:GenoArray
```
Return an array of all the non-missing alleles of a locus. Use
`miss = true` to include missing values.

----

### `unique_alleles`
```julia
unique_alleles(locus::T) where T<:GenotypeArray
```
Return an array of all the unique non-missing alleles of a locus.

----

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

----

### `copy`
```julia
copy(data::PopData)
```
Creates a copy of `PopData`.

----

### `drop_monomorphic`
    drop_monomorphic(data::PopData)
Return a `PopData` object omitting any monomorphic loci. Will inform you which loci were removed.


----

### `drop_monomorphic!`
    drop_monomorphic!(data::PopData)
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which loci were removed.

----

### `loci_dataframe`
    loci_dataframe(data::PopData)
Return a wide `DataFrame` of samples as columns, ommitting population information.

**Example**
```
julia> loci_dataframe(@nancycats)
9×237 DataFrame. Omitted printing of 232 columns
│ Row │ N215       │ N216       │ N217       │ N218       │ N219       │
│     │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │
├─────┼────────────┼────────────┼────────────┼────────────┼────────────┤
│ 1   │ missing    │ missing    │ (135, 143) │ (133, 135) │ (133, 135) │
│ 2   │ (136, 146) │ (146, 146) │ (136, 146) │ (138, 138) │ (140, 146) │
│ 3   │ (139, 139) │ (139, 145) │ (141, 141) │ (139, 141) │ (141, 145) │
│ 4   │ (116, 120) │ (120, 126) │ (116, 116) │ (116, 126) │ (126, 126) │
│ 5   │ (156, 156) │ (156, 156) │ (152, 156) │ (150, 150) │ (152, 152) │
│ 6   │ (142, 148) │ (142, 148) │ (142, 142) │ (142, 148) │ (142, 148) │
│ 7   │ (199, 199) │ (185, 199) │ (197, 197) │ (199, 199) │ (193, 199) │
│ 8   │ (113, 113) │ (113, 113) │ (113, 113) │ (91, 105)  │ (113, 113) │
│ 9   │ (208, 208) │ (208, 208) │ (210, 210) │ (208, 208) │ (208, 208) │
```

----

### `loci_matrix`
    loci_matrix(data::PopData)
Return a matrix of genotypes with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> loci_matrix(@nancycats)
237×9 Array{Union{Missing, Tuple{Int16,Int16}},2}:
 missing     (136, 146)  (139, 139)  …  (199, 199)  (113, 113)  (208, 208)
 missing     (146, 146)  (139, 145)     (185, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (141, 141)     (197, 197)  (113, 113)  (210, 210)
 (133, 135)  (138, 138)  (139, 141)     (199, 199)  (91, 105)   (208, 208)
 (133, 135)  (140, 146)  (141, 145)     (193, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (145, 149)  …  (193, 195)  (91, 113)   (208, 208)
 (135, 135)  (136, 146)  (139, 145)     (199, 199)  (105, 113)  (208, 208)
 (135, 143)  (136, 146)  (135, 149)     (193, 197)  (91, 91)    (208, 212)
 (137, 143)  (136, 146)  (139, 139)     (197, 197)  (105, 113)  (208, 212)
 (135, 135)  (132, 132)  (141, 145)     (197, 197)  (91, 105)   (208, 208)
 (137, 141)  (130, 136)  (137, 145)  …  (193, 199)  (91, 91)    (182, 182)
 (129, 133)  (130, 136)  (135, 145)     (193, 199)  (91, 113)   (182, 208)
 ⋮                                   ⋱                          
 (133, 135)  (136, 136)  (135, 139)  …  (199, 199)  (113, 113)  (182, 182)
 (133, 141)  (136, 136)  (135, 139)     (197, 197)  (113, 113)  (182, 208)
 (133, 141)  (130, 146)  (141, 141)     (191, 199)  missing     (208, 208)
 (123, 133)  (138, 138)  (141, 145)     (191, 197)  missing     (208, 208)
 (123, 133)  (138, 138)  (139, 139)     (197, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 139)  …  (197, 197)  missing     (208, 208)
 (133, 141)  (130, 136)  (139, 145)     (191, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 145)     (199, 199)  missing     (208, 220)
 (133, 143)  (130, 130)  (135, 145)     (197, 197)  missing     (208, 208)
 (135, 141)  (136, 144)  (143, 143)     (191, 197)  (113, 117)  (208, 208)
 (137, 143)  (130, 136)  (135, 145)  …  (193, 199)  (113, 117)  (208, 208)
 (135, 141)  (130, 146)  (135, 139)     (197, 197)  missing     (208, 208)
 ```



----

### `motivational_quote`
```julia
motivational_quote()
```
Randomly returns one of ~50 _somewhat_ motivational quotes.

----

### `nonmissing`
```julia
nonmissing(vec::T) where T<:AbstractArray
```
Convenience function to count the number of non-`missing` values in a vector.

    nonmissing(data::PopData, locus::String)
Convenience function to count the number of non-`missing` samples
at a locus.

----

### `nonmissings`
    nonmissings(vec1::AbstractVector, vec2::AbstractVector)
Return a vector of indices where neither input vectors have a `missing` value, i.e. an
intersection of the indices of their non-missing elements.

----

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

----

### `pairwise_pairs`
    pairwise_pairs(smp_names::Vector{String})
Given a vector of sample names, returns a vector of tuples of unique all x 
all combinations of sample pairs, excluding self-comparisons.

**Example**
```
julia> samps = ["red_1", "red_2", "blue_1", "blue_2"] ;
julia> pairwise_pairs(samps)
6-element Array{Tuple{String,String},1}:
 ("red_1", "red_2")
 ("red_1", "blue_1")
 ("red_1", "blue_2")
 ("red_2", "blue_1")
 ("red_2", "blue_2")
 ("blue_1", "blue_2")
```

----

### `phase`
    phase(data::PopData)
Return a `Vector` of length `ploidy` composed of allele matrices with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> mtx = phase(@nancycats)
2-element Array{Array{Union{Missing, Int16},2},1}:
 [missing 136 … 113 208; missing 146 … 113 208; … ; 137 130 … 113 208; 135 130 … missing 208]
 [missing 146 … 113 208; missing 146 … 113 208; … ; 143 136 … 117 208; 141 146 … missing 208]
julia> mtx[1]
237×9 Array{Union{Missing, Int16},2}:
    missing  136  139  116         156  142  199  113         208
    missing  146  139  120         156  142  185  113         208
 135         136  141  116         152  142  197  113         210
 133         138  139  116         150  142  199   91         208
 133         140  141  126         152  142  193  113         208
 135         136  145  120         150  148  193   91         208
 135         136  139  116         152  142  199  105         208
 135         136  135  120         154  142  193   91         208
 137         136  139  116         150  142  197  105         208
 135         132  141  120         150  148  197   91         208
 137         130  137  128         152  142  193   91         182
 129         130  135  126         144  140  193   91         182
   ⋮                                      ⋮                   
 133         136  135     missing  146  142  199  113         182
 133         136  135     missing  150  142  197  113         182
 133         130  141     missing  148  142  191     missing  208
 123         138  141     missing  148  142  191     missing  208
 123         138  139     missing  150  142  197     missing  208
 133         136  139     missing  150  142  197     missing  208
 133         130  139     missing  152  142  191     missing  208
 133         136  139     missing  150  142  199     missing  208
 133         130  135     missing  148  142  197     missing  208
 135         136  143     missing  144  142  191  113         208
 137         130  135     missing  150  142  193  113         208
 135         130  135     missing  150  142  197     missing  208
```

----

### `quickstart`
```julia
quickstart()
```
Prints helpful text of how to get started using PopGen.

----

### `reciprocal`
```julia
reciprocal(num::T) where T <: Signed
```
Returns the reciprocal (1/number) of a number. Will return `0` when 
the number is `0` instead of returning `Inf`.

----

### `strict_shuffle`
```julia
strict_shuffle(x::T) where T <: AbstractArray
```
Shuffle only the non-missing values of a Vector and return a copy of the vector,
keeping the `missing` values at their original locations.
Use `strict_shuffle!` to edit in-place instead of returning a copy.

----

### `strict_shuffle!`
```julia
strict_shuffle!(x::T)! where T <: AbstractArray
```
Shuffle only the non-missing values of a Vector, keeping the
`missing` values at their original locations. Use `strict_shuffle`
to return a copy instead of editing in-place.