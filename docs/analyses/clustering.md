---
id: clustering
title: Clustering
sidebar_label: Clustering
---

Usually the beginning of a study without prior population information requires guesstimating the number of clusters present in the data. This can be accomplished
using a number of methods, like K-means, K-mediods, Fuzzy-C Means, etc. PopGen.jl
extends several of the clustering algorithms available in [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) to work directly with PopData objects.

## The `cluster` wrapper
All of the clustering methods implemented in PopGen.jl (read below) can be accessed using a single function `cluster`.
```julia
cluster(::PopData, method::Function; matrixtype::Symbol, kwargs...)
```
A convenience wrapper to perform clustering on a `PopData` object determined by a designated `method`. The
chosen method must also be supplied with the appropriate keyword arguments for that method. For more information on 
a specific method, read more below or see its docstring in a Julia session with `?methodname` (e.g., `?kmediods`). The keyword argument `matrixtype` refers to which input matrix you would like to use for clustering, one of either `:pca` (default, principal components of the scaled allele frequencies) or `:freq` (scaled allele frequencies).

#### Clustering Methods
| Method Name | Method Type | Keyword Arguments |
|:---|:---|:---|
|`kmeans`| K-means++ | `k`, `iterations` |
|`kmedoids`| K-medoids | `k`, `iterations` |
|`hclust`| Hierarchical Clustering | `linkage`, `branchorder`, `distance` |
|`fuzzycmeans`| Fuzzy C-means | `c`, `fuzziness`, `iterations` |
|`dbscan`| Density-based Spatial Clustering of Applications with Noise (DBSCAN) | `k`, `iterations` |`radius`, `minpoints`, `distance` |

#### Examples
```julia
julia> cats = @nancycats;

julia> cluster(cats, kmeans, iterations = 100);

julia> cluster(cats, dbscan, matrixtype = :freq)
```

The results of these clustering methods can then be used for validation using any methods available in [Clustering.jl](https://github.com/JuliaStats/Clustering.jl).

:::info skipping the wrapper
Since the clustering methods are exported, you can technically skip the `cluster` wrapper and use any of the methods directly (e.g. `kmeans(PopData, k = 5)`), although `cluster()` is the preferred method.
:::

## Clustering Methods
### K-means
```julia
kmeans(data::PopData; k::Int64, iterations::Int64 = 100, matrixtype::Symbol = :pca)
```    
Perform Kmeans clustering (using Kmeans++ from [Arthur & Vassilvitskii 2007](http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf)) on a `PopData` object. Returns a `KmeansResult`
object. Use the keyword argument `iterations` (default: `100`) to set the maximum number of iterations allowed to
achieve convergence. Clustering is performed on the `matrixtype` principal components of the scaled allele frequencies (`:pca`),
or just the scaled allele frequencies themselves (`:freq`). In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `k`: the number of desired clusters, given as an `Integer`
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components of `:freq`
  - `:freq`: matrix of scaled allele frequencies 

**Example**
```julia
julia> cats = @nancycats ;

julia> km = kmeans(cats, k = 2)
```

----

### K-medoids
```julia
kmedoids(data::PopData; k::Int64, iterations::Int64 = 100, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
Perform K-medoids ([Kaufman & Rousseeuw, 1990]( https://doi.org/10.1002/9780470316801.ch2)) clustering on a `PopData` object. Returns a `KmedoidsResult`
object. Use the keyword argument `iterations` (default: `100`) to set the maximum number of iterations allowed to
achieve convergence. Clustering is performed on the `matrixtype` principal components of the scaled allele frequencies (`:pca`),
or just the scaled allele frequencies themselves (`:freq`). In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `k`: the number of desired clusters, given as an `Integer`
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. `sqeuclidean`, etc.)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components of `:freq`
  - `:freq`: matrix of scaled allele frequencies 

**Example**
```julia
julia> cats = @nancycats ;

julia> km = kmedoids(cats, k = 2, distance = sqeuclidean)
```

----

### Hierarchical Clustering
```julia
hclust(data::PopData; linkage::Symbol = :single, branchorder::Symbol = :r, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
Perform hierarchical clustering ([Bar-Joseph *et al.*, 2001](https://doi.org/10.1093/bioinformatics/17.suppl_1.S22)) on a PopData object. Returns an `Hclust` object, which contains many metrics but does not include cluster assignments. Use 
`cutree(::PopData, ::Hclust; krange...)` to compute the sample assignments for a range of `k` clusters. Clustering is performed on the `matrixtype` principal components of the scaled allele frequencies (`:pca`),
or just the scaled allele frequencies themselves (`:freq`). In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `linkage`: defines how the distances between the data points are aggregated into the distances between the clusters
  - `:single`: use the minimum distance between any of the cluster members (default)
  - `:average`: use the mean distance between any of the cluster members
  - `:complete`: use the maximum distance between any of the members
  - `:ward`: the distance is the increase of the average squared distance of a point to its cluster centroid after merging the two clusters
  - `:ward_presquared`: same as `:ward`, but assumes that the distances in the distance matrix are already squared.
- `branchorder`: algorithm to order leaves and branches (default: `:r`)
  - `:r`: ordering based on the node heights and the original elements order (compatible with R's hclust)
  - `:optimal`: branches are ordered to reduce the distance between neighboring leaves from separate branches using the "fast optimal leaf ordering" [algorithm](https://doi.org/10.1093/bioinformatics/17.suppl_1.S22)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. `sqeuclidean`, etc.)
- `matrixtype`: type of input matrix (default: `:pca`)
  - `:pca`: matrix of Principal Components of `:freq`
  - `:freq`: matrix of allele frequencies

#### cutree
```julia
cutree(::PopData, hcres::Hclust; krange::UnitRange{Int64}, height::Union{Int64, Nothing} = nothing)
cutree(::PopData, hcres::Hclust; krange::Vector{Int64}, height::Union{Int64, Nothing} = nothing)
```
An expansion to the `Clustering.cutree` method (from Clustering.jl) that performs cluster assignments over `krange`
on the `Hclust` output from `hclust()`. Returns a `DataFrame` of sample names and columns corresponding to assignments 
per k in `krange`. The `PopData` object is used only for retrieving the sample names.

**Keyword Arguments**
- `krange`: the number of desired clusters, given as a vector (ex. `[2,4,5]`) or range (`2:5`)
- `h::Integer`: the height at which the tree is cut (optional) 

**Example**
```julia
julia> cats = @nancycats ;

julia> hca = hclust(cats, branchorder = :optimal) ;

julia> cutree(cats, hca, krange = 2:5)
```

----

### Fuzzy C-means
```julia
fuzzycmeans(data::PopData; c::Int64, fuzziness::Int64 = 2, iterations::Int64 = 100, matrixtype::Symbol = :pca)
```
Perform Fuzzy C-means clustering ([Bezdek *et al.* 1984](https://doi.org/10.1016/0098-3004(84)90020-7)) on a PopData object. Returns a `FuzzyCMeansResult` object, which contains the assignment weights in the `.weights` field. Clustering is performed on the `matrixtype` principal components of the scaled allele frequencies (`:pca`),
or just the scaled allele frequencies themselves (`:freq`). In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `c`: the number of desired clusters, given as an `Integer`
- `fuzziness::Integer`: clusters' fuzziness, must be >1 (default: `2`)
  - a fuzziness of 2 is common for systems with unknown numbers of clusters
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components of `:freq`
  - `:freq`: matrix of scaled allele frequencies 

**Example**
```julia
julia> cats = @nancycats ;

julia> fuzzycats = fuzzycmeans(cats, c = 5) ;
```

----

### DBSCAN
```julia
dbscan(::PopData; radius::Float64, minpoints::Int64 = 2, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
  Perform Density-based Spatial Clustering of Applications with Noise (DBSCAN: [Ester *et al.* 1996](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.121.9220))
on a PopData object. Returns a `DbscanResult` object, which contains the assignments in the
`.assignments` field. Clustering is performed on the `matrixtype` principal components of the scaled allele frequencies (`:pca`),
or just the scaled allele frequencies themselves (`:freq`). In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `radius::Float64`: the radius of a point neighborhood
- `minpoints::Int`: the minimum number of a core point neighbors (default: `2`)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. `sqeuclidean`, etc.)
- `matrixtype`: type of input matrix (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of allele frequencies 

**Example**
```julia
julia> cats = @nancycats ;

julia> fuzzycats = dbscan(cats, radius = 0.5) ;
```

----------
## Acknowledgments
Much of the heavy lifting within these clustering methods are actually the result of the
amazing authors and contributors of [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) and the Principal Component Analysis available from [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl).