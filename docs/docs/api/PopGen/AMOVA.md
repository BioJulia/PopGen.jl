---
id: amova
title: AMOVA.jl
sidebar_label: AMOVA.jl
---

## PopGen.jl/src/AMOVA.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 📦 allelicdistance
```julia
allelicdistance(x, y)
```

### 📦 _missinglocusfilter
```julia
_missinglocusfilter(data::PopData, cutoff::Float64)
```

### 🔵 AMOVA
```julia
AMOVA(data::PopData; kwargs...) = amova(data, kwargs...)
```
This function is needed for `pairwisefst` dispatch, so it doesn't use the full AMOVA calculation.

### 🔵 amova
```julia
amova(data::PopData; hierarchy::String, missing_cutoff::Union{Nothing,Float64} = 0.05)
```

**Arguments**
- `data`: A `PopData` object

**Keyword Arguments**

- `hierarchy`: A `String` of what metadata column has the grouping type you would like to use. This is likely the `popultion` column in `sampleinfo`, since heirarchical clustering is not yet supported. The `heirarchy` API is very likely to change once that feature is implemented.
- `missing_cutoff`: A threshold to use for filtering out missing data. The default is `0.05`, meaning remove loci with greater than 5% missing data. It accepts `nothing` to indicate you want all data kept.

The result is of type `AMOVAResult`, which allows it to be nicely printed as a familiar AMOVA/ANOVA table. You can easily access the fields of the result using dot-indexing, _e.g._ `result.SS`. Note that sigma-squared uses fancy characters, which will require Julia's built-in ASCII character rendering: `\sigma<tab>\^2<tab>`, where `<tab>` is the tab key on your keyboard.

### 🔵 AMOVAResult
```julia
struct AMOVAResult
    source::Vector{String}
    df::Vector{Int}
    SS::Vector{Float64}
    MS::Vector{Float64}
    σ²::Vector{Float64}
    FST::Float64
end
```
### 🔵 show
```julia
Base.show(io::IO, result::AMOVAResult)
```
Borrowed heavily from [StatsModels.jl](https://github.com/JuliaStats/StatsModels.jl/blob/6f19ecca344b4dc99c41d0e70976f306a3dd9e72/src/lrtest.jl#L130) to have consistent style with expected
model outputs. 
