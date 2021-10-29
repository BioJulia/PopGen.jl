---
id: utils
title: Utils.jl
sidebar_label: Utils.jl
---

## PopGen.jl/src/Utils.jl
❗ => not exported | 
🔵 => exported by PopGen.jl

### ❗_adjacency_matrix
```julia
_adjacency_matrix(data::PopData)
```

### ❗_p_adjust
```julia
_p_adjust(pvals::Vector{T}, method::String) where T <: Union{Missing, <:AbstractFloat}
```
Modification to `MultipleTesting.adjust` to include `missing` values in the
returned array. See MultipleTesting.jl docs for full more detailed information.

**Example**
```
julia> _p_adjust([0.1, 0.01, 0.005, 0.3], "bh")
```
**`correction` methods (case insensitive)**
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
"""


### ❗feature_req
```julia
feature_req()
```
Returns the text: `Please open an Issue or Pull Request on https://www.github.com/biojulia/PopGen.jl if you would like this feature implemented`