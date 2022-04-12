#struct ClusteringResults
#    method::String
#    assignments::DataFrame
#    costs::DataFrame
#    other::DataFrame
#    centers::NamedTuple
#end

#function Base.show(io::IO, data::ClusteringResults)
#    printstyled(io, "$(data.method) Clustering Results\n", bold = true)
#    println(io, "  K values:          " * join(string.(data.other.k), ", "))
#    println(io, "  Iterations per K:  " * join(string.(data.other.iterations), ", "))
#    println(io, "  Convergence per K: " * join(replace(x -> x == "true" ? "T" : "F", string.(data.other.converged)), ", "))
#    println(io, "  Total Cost per K:  " * join(string.(round.(data.other.totalcost, digits = 3)), ", "))
#    println(io, "\n  Available fields to inspect: assignments, costs, centers, other")
#end


"""
    kmeans(data::PopData; k::Int64, iterations::Int64 = 100, matrixtype::Symbol = :pca)
    
Perform Kmeans clustering (using Kmeans++) on a `PopData` object. Returns a `ClusteringResults`
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
"""
function kmeans(data::PopData; k::Int64, iterations::Int64 = 100, matrixtype::Symbol = :pca)
    mtx = matrixtype == :pca ?
        pca(data, center = false, scale = true).proj |> permutedims : matrixtype == :freq ?
        _allelematrix(data, center = false, scale = true) :
        throw(ArgumentError("matrixtype :$matrixtype invalid, choose between :pca or :freq"))
    kmeans(mtx, k, maxiter = iterations)
    #idx = 1:length(out)
    #assn = DataFrame([getproperty(out[i], :assignments) for i in idx], Symbol.(krange))
    #insertcols!(assn, 1, :name => unique(data.genodata.name))
    #centers = NamedTuple{Tuple(Symbol.(krange))}(Tuple(getproperty(out[i], :centers) for i in idx))
    #costs = DataFrame([getproperty(out[i], :costs) for i in idx], Symbol.(krange))
    #insertcols!(costs, 1, :name => assn.name)
    #other = 
    #    DataFrame(
    #        :k => krange,
    #        :iterations => [getproperty(out[i], :iterations) for i in idx],
    #        :converged => [getproperty(out[i], :converged) for i in idx],
    #        :totalcost => [getproperty(out[i], :totalcost) for i in idx],
    #        :counts => [getproperty(out[i], :counts) for i in idx],
    #        :wcounts => [getproperty(out[i], :wcounts) for i in idx],
    #        :cweights => [getproperty(out[i], :cweights) for i in idx]
    #    )
    #ClusteringResults("K-means++", assn, costs, other, centers)
end


"""
    kmedoids(data::PopData; krange::Int64, iterations::Int64 = 100, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)

Perform Kmedoids clustering on a `PopData` object. Returns a `ClusteringResults`
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
"""
function kmedoids(data::PopData; k::Int64, iterations::Int64 = 100, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
    mtx = matrixtype == :pca ?
        pairwise(distance, pca(data, center = false, scale = true).proj |> permutedims, dims = 2) : matrixtype == :freq ?
        pairwise(distance, _allelematrix(data, center = false, scale = true), dims = 1) : 
        throw(ArgumentError("matrixtype :$matrixtype invalid, choose between :pca or :freq"))
    kmedoids(mtx, k, maxiter = iterations)
end


"""
    hclust(data::PopData; linkage::Symbol = :single, branchorder::Symbol = :r, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)

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
"""
function hclust(data::PopData; linkage::Symbol = :single, branchorder::Symbol = :r, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
    mtx = matrixtype == :pca ?
        pairwise(distance, pca(data, center = false, scale = true).proj |> permutedims, dims = 2) : matrixtype == :freq ?
        pairwise(distance, _allelematrix(data, center = false, scale = true), dims = 1) : 
        throw(ArgumentError("matrixtype :$matrixtype invalid, choose between :pca or :freq"))
    hclust(mtx, linkage = linkage, branchorder = branchorder)
end

"""
    cutree(::PopData, hcres::Hclust; krange::UnitRange{Int64}, height::Union{Int64, Nothing} = nothing)
    cutree(::PopData, hcres::Hclust; krange::Vector{Int64}, height::Union{Int64, Nothing} = nothing)

An expansion to the `Clustering.cutree` method (from Clustering.jl) that performs cluster assignments over `krange`
on the `Hclust` output from `hclust()`. Returns a `DataFrame` of sample names and columns corresponding to assignments 
per k in `krange`. The `PopData` object is used only for retrieving the sample names.

**Keyword Arguments**
- `krange`: the number of desired clusters, given as a vector (ex. `[2,4,5]`) or range (`2:5`)
- `h::Integer`: the height at which the tree is cut (optional) 

"""
function cutree(data::PopData, hcres::Hclust; krange::Union{UnitRange{Int64},Vector{Int64}}, height::Union{Int64, Nothing} = nothing)
    tmp = map(i -> cutree(hcres, k = i, h = height), krange)
    out = DataFrame([tmp[i] for i in 1:length(krange)], Symbol.(krange))
    insertcols!(out, 1, :name => unique(data.genodata.name))
    return out
end

"""
    fuzzycmeans(data::PopData; c::Int64, fuzziness::Int64 = 2, iterations::Int64 = 100, matrixtype::Symbol = :pca)
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
"""
function fuzzycmeans(data::PopData; c::Int64, fuzziness::Int64 = 2, iterations::Int64 = 100, matrixtype::Symbol = :pca)
    mtx = matrixtype == :pca ?
        pca(data, center = false, scale = true).proj |> permutedims : matrixtype == :freq ?
        _allelematrix(data, center = false, scale = true) :
        throw(ArgumentError("matrixtype :$matrixtype invalid, choose between :pca or :freq"))
    fuzzy_cmeans(mtx, c, fuzziness, maxiter = iterations)
end

"""
    dbscan(data::PopData; radius::Float64, minpoints::Int64 = 2, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)

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
"""
function dbscan(data::PopData; radius::Float64, minpoints::Int64 = 2, distance::PreMetric = euclidean, matrixtype::Symbol = :pca)
    mtx = matrixtype == :pca ?
        pairwise(distance, pca(data, center = false, scale = true).proj |> permutedims, dims = 2) : matrixtype == :freq ?
        pairwise(distance, _allelematrix(data, center = false, scale = true), dims = 1) : 
        throw(ArgumentError("matrixtype :$matrixtype invalid, choose between :pca or :freq"))
    dbscan(mtx, radius, minpoints)
end