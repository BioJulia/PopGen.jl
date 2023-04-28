---
id: clustering
title: Clustering.jl
sidebar_label: Clustering.jl
---

## PopGen.jl/src/Clustering.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|


### 🔵 cutree
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

----
### 🔵 cluster
```julia
cluster(::PopData, method::Function ; kwargs)
```
A convenience wrapper to perform clustering on a `PopData` object determined by a designated `method` (see below). The
chosen method must also be supplied with the appropriate keyword arguments for that method. For more information on 
a specific method, see its docstring with `?methodname`

**Clustering Methods**
- `kmeans`: K-means++ clustering
  - kwargs: `k`, `iterations`, `matrixtype`
- `kmedoids`: K-medoids clustering
  - kwargs: `k`, `iterations`, `distance`, `matrixtype`
- `hclust`: Hierarchical clustering
  - kwargs: `linkage`, `branchorder`, `distance`, `matrixtype`
- `fuzzycmeans`: Fuzzy C-means lustering
  - kwargs: `c`, `fuzziness`, `iterations`, `matrixtype`
- `dbscan`: Density-based Spatial Clustering of Applications with Noise (DBSCAN)
  - kwargs: `radius`, `minpoints`, `distance`, `matrixtype`

----

### 🔵 dbscan
```julia
dbscan(::PopData; radius::Float64, minpoints::Int64 = 2, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
An expansion of `Clustering.dbscan` (from Clustering.jl) to perform Density-based Spatial Clustering of Applications with Noise (DBSCAN)
on a PopData object. This is a convenience method which converts the `PopData` object to either an allele frequency or PCA matrix, and performs
DBSCAN clustering on the distance matrix of that. Returns a `DbscanResult` object, which contains the assignments in the
`.assignments` field.

**Keyword Arguments**
- `radius::Float64`: the radius of a point neighborhood
- `minpoints::Int`: the minimum number of a core point neighbors (default: `2`)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. sqeuclidean, etc.)
- `matrixtype`: type of input matrix (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of allele frequencies 

----

### 🔵 fuzzycmeans
```julia
fuzzycmeans(data::PopData; c::Int64, fuzziness::Int64 = 2, iterations::Int64 = 100, matrixtype::Symbol = :pca)
```
An expansion of `Clustering.fuzzy_cmeans` (from Clustering.jl) to perform Fuzzy C-means clustering on a PopData object. This is a convenience method
which converts the `PopData` object to either an allele frequency or PCA matrix, and performs Fuzzy C-means
clustering on the Euclidean distance matrix of that. Returns a `FuzzyCMeansResult` object, which contains the assignment weights in the
`.weights` field.

**Keyword Arguments**
- `c`: the number of desired clusters, given as an `Integer`
- `fuzziness::Integer`: clusters' fuzziness, must be >1 (default: `2`)
  - a fuzziness of 2 is common for systems with unknown numbers of clusters
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of scaled allele frequencies 

----

### 🔵 hclust
```julia
hclust(data::PopData; linkage::Symbol = :single, branchorder::Symbol = :r, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
An expansion of `Clustering.hclust` (from Clustering.jl) to perform hierarchical clustering on a PopData object. This is a convenience method
which converts the `PopData` object to either an allele frequency or PCA matrix, converts that into a distance matrix, and performs hierarchical
clustering on that distance matrix. Returns an `Hclust` object, which contains many metrics but does not include cluster assignments. Use 
`cutree(::PopData, ::Hclust; krange...)` to compute the sample assignments for a range of `k` clusters.

**Keyword Arguments**
- `linkage`: defines how the distances between the data points are aggregated into the distances between the clusters (default: `:single`)
  - `:single`: use the minimum distance between any of the cluster members
  - `:average`: use the mean distance between any of the cluster members
  - `:complete`: use the maximum distance between any of the members
  - `:ward`: the distance is the increase of the average squared distance of a point to its cluster centroid after merging the two clusters
  - `:ward_presquared`: same as `:ward`, but assumes that the distances in the distance matrix are already squared.
- `branchorder`: algorithm to order leaves and branches (default: `:r`)
  - `:r`: ordering based on the node heights and the original elements order (compatible with R's hclust)
  - `:optimal`: branches are ordered to reduce the distance between neighboring leaves from separate branches using the "fast optimal leaf ordering" [algorithm](https://doi.org/10.1093/bioinformatics/17.suppl_1.S22)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. sqeuclidean, etc.)
- `matrixtype`: type of input matrix (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of allele frequencies 

----

### 🔵 kmeans
```julia
kmeans(data::PopData; k::Int64, iterations::Int64 = 100, matrixtype::Symbol = :pca)
```    
Perform Kmeans clustering (using Kmeans++) on a `PopData` object. Returns a `KmeansResult`
object. Use the keyword argument `iterations` (default: 100) to set the maximum number of iterations allowed to
achieve convergence. Interally, kmeans clustering is performed on either the principal components of the scaled allele frequencies,
or just the scaled allele frequencies themselves. In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `k`: the number of desired clusters, given as an `Integer`
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of scaled allele frequencies 

**Example**
```julia
julia> cats = @nancycats ;

julia> km = kmeans(cats, k = 2)
```

----

### 🔵 kmedoids
```julia
kmedoids(data::PopData; k::Int64, iterations::Int64 = 100, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
```
Perform Kmedoids clustering on a `PopData` object. Returns a `KmedoidsResult`
object. Use the keyword argument `iterations` (default: 100) to set the maximum number of iterations allowed to
achieve convergence. Interally, kmeans clustering is performed on either the principal components of the scaled allele frequencies,
or just the scaled allele frequencies themselves. In both cases, `missing` values are replaced by the global mean allele frequency.

**Keyword Arguments**
- `k`: the number of desired clusters, given as an `Integer`
- `iterations::Int64`: the maximum number of iterations to attempt to reach convergence (default: `100`)
- `distance`: type of distance matrix to calculate on `matrixtype` (default: `euclidean`)
  - see [Distances.jl](https://github.com/JuliaStats/Distances.jl) for a list of options (e.g. sqeuclidean, etc.)
- `matrixtype`: type of input matrix to compute (default: `:pca`)
  - `:pca`: matrix of Principal Components
  - `:freq`: matrix of scaled allele frequencies 

----

### 🔵 show
```julia
Base.show(io::IO, data::KMeansResults)
```
