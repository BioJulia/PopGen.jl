---
id: kmeans
title: K-Means Clustering
sidebar_label: K-Means Clustering
---

Usually the beginning of a study without prior population information requires guesstimating the number of clusters present in the data.
One way to accomplish this is with K-means clustering, an unsupervised clustering algorithm that clusters data based on similarity. The
PopGen.jl implementation of K-means clustering uses the K-means ++ algorithm ([Arthur & Vassilvitskii 2007](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf)) under the hood, as implemented in `Clustering.jl` ([link](https://github.com/JuliaStats/Clustering.jl)).


### `kmeans`
```julia
# krange as a vector : [1,3,6,7]
kmeans(data::PopData; krange::Vector{Int64}, iterations::Int64 = 100)
# krange as a range: 2:9
kmeans(data::PopData; krange::UnitRange{Int64}, iterations::Int64 = 100)
# without keywords
kmeans(data::PopData, krange::Union{UnitRange{Int64},Vector{Int64}}, iterations::Int64 = 100)
```

Perform PCA-based Kmeans clustering (using Kmeans++) on a `PopData` object. Returns a `KMeansResults`
object storing the results of each `kmeans()` run across a range of K values (`krange`, default: `2:nsamples/10`).
Use the keyword argument `iterations` (default: 100) to set the maximum number of iterations allowed to
achieve convergence. Interally, kmeans clustering is performed on the principal components of the scaled allele frequencies 
matrix, where `missing` values are replaced by the global mean allele frequency.

**Example**
```julia
julia> cats = @nancycats ;
julia> km = kmeans(cats, krange = 2:7)
K-Means(++) Clustering Results
  K values:          2, 3, 4, 5, 6, 7
  Iterations per K:  2, 4, 6, 5, 5, 7
  Convergence per K: T, T, T, T, T, T
  Total Cost per K:  91.41, 90.679, 89.54, 89.287, 88.284, 87.957
  Available fields to inspect: assignments, costs, centers, other
  
julia> km.assignments
237×7 DataFrame
 Row │ name     2      3      4      5      6      7     
     │ String   Int64  Int64  Int64  Int64  Int64  Int64
─────┼───────────────────────────────────────────────────
   1 │ N215         1      1      1      2      3      3
   2 │ N216         1      2      1      2      3      6
   3 │ N217         1      3      1      1      3      3
   4 │ N218         1      3      1      1      3      6
   5 │ N219         1      2      1      2      3      3
   6 │ N220         1      2      1      2      3      1
  ⋮  │    ⋮       ⋮      ⋮      ⋮      ⋮      ⋮      ⋮
 233 │ N296         1      1      1      2      3      6
 234 │ N297         1      1      3      2      4      3
 235 │ N281         1      2      1      2      3      3
 236 │ N289         1      3      3      2      3      4
 237 │ N290         1      1      1      2      3      6
```     

