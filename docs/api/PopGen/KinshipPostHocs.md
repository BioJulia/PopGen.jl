---
id: kinshipposthocs
title: KinshipPostHocs.jl
sidebar_label: KinshipPostHocs.jl
---
## PopGen.jl/src/Kinship/KinshipPostHocs.jl
📦  => not exported | 
🔵 => exported by PopGen.jl

### 📦 sig_within
```julia
sig_within(data::PopData, results::DataFrame, population::String, iterations::Int = 20000)
```

### 🔵 kinshipposthoc
```julia
kinshipposthoc(::PopData, results::Union{DataFrame, NamedTuple}; iterations::Int)
```
Performs a posthoc analysis using the resulting DataFrame or NamedTuple
from `kinship()`. This analysis uses permutations to test if a population has
significantly higher within-population kinship than between-population kinship.
The `results` object must have been generated from the provided `PopData`. Use `iterations = `
to specify the number of iterations for the permutation tests (default = `20000`). **Recommended**
that you use `MultipleTesting.jl` to correct resulting P-values.

**Example**
```julia
julia> cats = @nancycats ;
julia> rel_out = kinship(cats, method = [Ritland, Moran], iterations = 100);
julia> kinshipposthoc(cats, rel_out)
17x3 DataFrame
 Row │ population  Ritland_P  Moran_P
     │ String      Float64    Float64
─────┼────────────────────────────────
   1 │ 1              5.0e-5   5.0e-5
   2 │ 2              5.0e-5   5.0e-5
   3 │ 3              5.0e-5   5.0e-5
   4 │ 4              5.0e-5   5.0e-5
   5 │ 5              5.0e-5   5.0e-5
   6 │ 6              5.0e-5   5.0e-5
   7 │ 7              5.0e-5   5.0e-5
   8 │ 8              5.0e-5   5.0e-5
   9 │ 9              5.0e-5   5.0e-5
  10 │ 10             5.0e-5   5.0e-5
  11 │ 11             5.0e-5   5.0e-5
  12 │ 12             5.0e-5   5.0e-5
  13 │ 13             5.0e-5   5.0e-5
  14 │ 14             5.0e-5   5.0e-5
  15 │ 15             5.0e-5   5.0e-5
  16 │ 16             5.0e-5   5.0e-5
  17 │ 17             5.0e-5   5.0e-5
```