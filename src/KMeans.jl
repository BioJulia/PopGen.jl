struct KMeansResults
    assignments::DataFrame
    costs::DataFrame
    other::DataFrame
    centers::NamedTuple
end

function Base.show(io::IO, data::KMeansResults)
    printstyled(io, "K-Means(++) Clustering Results\n", bold = true)
    println(io, "  K values:          " * join(string.(data.other.k), ", "))
    println(io, "  Iterations per K:  " * join(string.(data.other.iterations), ", "))
    println(io, "  Convergence per K: " * join(replace(x -> x == "true" ? "T" : "F", string.(data.other.converged)), ", "))
    println(io, "  Total Cost per K:  " * join(string.(round.(data.other.totalcost, digits = 3)), ", "))
    println(io, "\n  Available fields to inspect: assignments, costs, centers, other")
end


"""
    kmeans(data::PopData; krange::Vector{Int64}, iterations::Int64 = 100)
    kmeans(data::PopData; krange::UnitRange{Int64}, iterations::Int64 = 100)
    kmeans(data::PopData, krange::Union{UnitRange{Int64},Vector{Int64}}, iterations::Int64 = 100)

Perform PCA-based Kmeans clustering (using Kmeans++) on a `PopData` object. Returns a `KMeansResults`
object storing the results of each `Clustering.kmeans()` run across a range of K values (`krange`, default: 2:nsamples/10).
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
     │ String7  Int64  Int64  Int64  Int64  Int64  Int64
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
                                         226 rows omitted
```
"""
function kmeans(data::PopData; krange::Union{UnitRange{Int64},Vector{Int64}} = 2:(size(data)[1] ÷ 10), iterations::Int64 = 100)
    kmeans(data, krange, iterations)
end

function kmeans(data::PopData, krange::Union{UnitRange{Int64},Vector{Int64}} = 2:(size(data)[1] ÷ 10), iterations::Int64 = 100)
    pc = pca(data, center = false, scale = true).proj |> permutedims
    out = map(i -> kmeans(pc, i, maxiter = iterations), krange)
    idx = 1:length(out)
    assn = DataFrame([getproperty(out[i], :assignments) for i in idx], Symbol.(krange))
    insertcols!(assn, 1, :name => unique(data.genodata.name))
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
