---
id: pairwisefst
title: PairwiseFST.jl
sidebar_label: PairwiseFST.jl
---

## PopGen.jl/src/FStatistics/PairwistFST.jl
❗ => not exported | 
🔵 => exported by PopGen.jl

### 🔵 Base.show
```julia
Base.show(io::IO, data::PairwiseFST)
```

-----

### 🔵 pairwisefst
```julia
pairwise_fst(data::PopData; method::String, iterations::Int64)
```
Calculate pairwise FST between populations in a `PopData` object. Set `iterations` 
to a value greater than `0` to perform a single-tailed permutation test to obtain
P-values of statistical significance.
#### Methods:
- `"Nei87"`: Nei (1987) method
- `"WC84"` : Weir-Cockerham (1984) method (default)
**Examples**

```julia
data = @nancycats
wc = pairwise_fst(data, method = "WC84")
wc_sig = pairwise_fst(data, iterations = 1000)
```
