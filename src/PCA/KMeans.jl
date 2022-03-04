struct KMeansResults
    assignments::DataFrame
    costs::DataFrame
    other::DataFrame
    centers::NamedTuple
end

function Base.show(io::IO, data::KMeansResults)
    println(io, "  KMeans results for K = [" * join(string.(data.other.k), ", ") * "]")
    println(io, "  Iterations per K: " * join(string.(data.other.iterations), ", "))
    println(io, "  Convergence per K: " * join(string.(data.other.converged), ", "))
    println(io, "  Total Cost per K: " * join(string.(round.(data.other.totalcost, digits = 3)), ", "))
    println(io, "\n  Available fields to inspect: assignments, costs, centers, other")
end


"""
    kmeans(data::PopData; krange::Vector{Int64}, iterations::Int64 = 100)
    kmeans(data::PopData; krange::UnitRange{Int64}, iterations::Int64 = 100)

Perform PCA-based Kmeans clustering (using Kmeans++) on a `PopData` object. Returns a `KMeansResults`
object storing the results of each `Clustering.kmeans()` run across a range of K values (`krange`, default: 2:nsamples/10).
Use the keyword argument `iterations` (default: 100) to set the maximum number of iterations allowed to
achieve convergence. Interally, kmeans clustering is performed on the principal components of the scaled allele frequencies matrix, where `missing` values are replaced by the global mean allele frequency.

**Example**
```julia
julia> cats = @nancycats ;

julia> km = kmeans(cats, krange = 2:7)
KMeans results for K = [2, 3, 4, 5, 6, 7]
    Iterations per K: 2, 3, 5, 6, 6, 5
    Convergence per K: true, true, true, true, true, true
    Total Cost per K: 91.084, 90.768, 89.722, 88.739, 88.295, 87.553
  
    Available fields to inspect: assignments, costs, centers, other

julia> km.assignments
237×6 DataFrame
 Row │ 2      3      4      5      6      7     
     │ Int64  Int64  Int64  Int64  Int64  Int64 
─────┼──────────────────────────────────────────
   1 │     2      3      3      3      3      4
   2 │     2      3      3      3      3      4
   3 │     2      2      3      3      1      6
   4 │     2      3      1      3      3      4
   5 │     2      3      3      3      3      4
   6 │     2      3      3      3      3      6
   7 │     2      3      3      3      3      4
   8 │     2      3      3      5      3      6
   9 │     2      3      3      5      6      4
  10 │     2      3      3      3      1      4
  ⋮  │   ⋮      ⋮      ⋮      ⋮      ⋮      ⋮
 229 │     2      3      1      3      1      4
 230 │     2      3      1      3      3      4
 231 │     2      3      3      3      3      4
 232 │     2      3      3      3      1      4
 233 │     2      3      3      3      3      4
 234 │     2      3      3      3      6      4
 235 │     2      3      3      3      1      2
 236 │     2      2      3      3      3      4
 237 │     2      3      3      3      3      4
                                218 rows omitted
```
"""
function Clustering.kmeans(data::PopData; krange::Union{UnitRange{Int64},Vector{Int64}} = 2:(size(data)[1] ÷ 10), iterations::Int64 = 100)
    pc = pca(data, center = false, scale = true).proj |> permutedims
    out = map(i -> Clustering.kmeans(pc, i, maxiter = iterations), krange)
    idx = 1:length(out)
    assn = DataFrame([getproperty(out[i], :assignments) for i in idx], Symbol.(krange))
    centers = NamedTuple{Tuple(Symbol.(krange))}(Tuple(getproperty(out[i], :centers) for i in idx))
    costs = DataFrame([getproperty(out[i], :costs) for i in idx], Symbol.(krange))
    other = 
        DataFrame(
            :k => krange,
            :iterations => [getproperty(out[i], :iterations) for i in idx],
            :converged => [getproperty(out[i], :converged) for i in idx],
            :totalcost => [getproperty(out[i], :totalcost) for i in idx],
            :counts => [getproperty(out[i], :counts) for i in idx],
            :wcounts => [getproperty(out[i], :wcounts) for i in idx],
            :cweights => [getproperty(out[i], :cweights) for i in idx]
        )
    KMeansResults(assn, costs, other, centers)
end
