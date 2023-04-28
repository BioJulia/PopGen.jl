---
id: amova
title: Analysis of Molecular Variance (AMOVA)
sidebar_label: AMOVA
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />

The Analysis of Molecular Variance, or AMOVA, was originally described in 
[Excoffier et _al._ 1992](https://doi.org/10.1093/genetics/131.2.479) and
is now a common method in the population genetics toolkit. The method 
compares within-population and across-population variance in allele 
frequencies to derive a metric of genetic differentiation. It has also been 
extended to calculate pairwise $F_{ST}$, which is the default method for 
`pairwisefst()` since PopGen.jl v0.10. What makes this method great is that
it can explore hierarchical population structure. That is, you can see if
there is significant structure as a result of (_e.g._) localities nested in
regions, or localities nested in regions nested in years, _etc_. Sadly, hierarchical structure is **not yet supported**, because I need time
to figure out the proper math for it. Help is always appreciated!

## AMOVA
You can perform an AMOVA using the `amova` function:
```julia
amova(data::PopData; hierarchy::String, missing_cutoff::Union{Nothing,Float64} = 0.05)
```
### Arguments
- `data`: A `PopData` object

### Keyword Arguments
- `hierarchy`: A `String` of what metadata column has the grouping type you would like to use. This is likely the `popultion` column in `sampleinfo`, since heirarchical clustering is not yet supported. The `heirarchy` API is very likely to change once that feature is implemented.
- `missing_cutoff`: A threshold to use for filtering out missing data. The default is `0.05`, meaning remove loci with greater than 5% missing data. It accepts `nothing` to indicate you want all data kept.

The result is of type `AMOVAResult`, which allows it to be nicely printed as a familiar AMOVA/ANOVA table. You can easily access the fields of the result using dot-indexing, _e.g._ `result.SS`. Note that sigma-squared uses fancy characters, which will require Julia's built-in ASCII character rendering: `\sigma<tab>\^2<tab>`, where `<tab>` is the tab key on your keyboard.
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


#### example
```julia
julia> cats = @nancycats;

julia> julia> amova(cats, hierarchy = "population")
[ Info: removing 2 loci with >5.0% missing data
Analysis of Molecular Variance
─────────────────────────────────────────────────────────
            Source   DF        SS      MS      σ²     FST
─────────────────────────────────────────────────────────
             Total  236  753.7099  8.8930  0.4409  0.1369
 Among populations   16       142  2.7792  2.7792        
Within populations  220       611                        
─────────────────────────────────────────────────────────
```

