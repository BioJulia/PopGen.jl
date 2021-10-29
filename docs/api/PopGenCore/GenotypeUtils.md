---
id: genotypeutils
title: GenotypeUtils.jl
sidebar_label: GenotypeUtils.jl
---
## PopGenCore.jl/src/Utils/GenotypeUtils.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl


### 🟪 alleles
```julia
allelecount(locus::T) where T<:GenoArray
```
Return the number of unique alleles present at a locus.

----
### 🟪 alleles
```julia
alleles(locus::T) where T<:GenoArray
```
Return an array of all the non-missing alleles of a locus.

----
### 🟪 alleles
```julia
alleles(locus::T, miss::Bool = false) where T<:GenoArray
```
Return an array of all the non-missing alleles of a locus. Use the second positional
argument as `true` to include missing values.

----
### 🟪 uniquealleles
```julia
uniquealleles(locus::T) where T<:GenoArray
```
Return an array of all the unique non-missing alleles of a locus.

----
### 🟪 locidataframe
```julia
locidataframe(data::PopData)
```
Return a wide `DataFrame` of samples as columns, ommitting population information.

**Example**
```
julia> locidataframe(@nancycats)
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
### 🟪 locimatrix
```julia
locimatrix(data::PopData)
```
Return a matrix of genotypes with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples.

**Example**
```
julia> locimatrix(@nancycats)
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
### 🟪 phasedmatrix
```julia
phasedmatrix(data::PopData)
```
Return a `Vector` of length `ploidy` composed of allele matrices with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> mtx = phasedmatrix(@nancycats)
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