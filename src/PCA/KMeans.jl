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
    println(io, "\n  Available fields to inspect: assignments, centers, costs, other")
end

function Clustering.kmeans(data::PopData; krange::Union{UnitRange{Int64}, Vector{Int}} = 0:0, iterations::Int64 = 100)
    kr = krange == 0:0 ? (2:(size(data)[1] ÷ 10)) : krange
    pc = pca(data, center = false, scale = true).proj |> permutedims
    out = map(i -> Clustering.kmeans(pc, i, maxiter = iterations), kr)
    idx = 1:length(out)
    assn = DataFrame([getproperty(out[i], :assignments) for i in idx], Symbol.(kr))
    centers = NamedTuple{Tuple(Symbol.(kr))}(Tuple(getproperty(out[i], :centers) for i in idx))
    costs = DataFrame([getproperty(out[i], :costs) for i in idx], Symbol.(kr))
    other = 
        DataFrame(
            :k => kr,
            :iterations => [getproperty(out[i], :iterations) for i in idx],
            :converged => [getproperty(out[i], :converged) for i in idx],
            :totalcost => [getproperty(out[i], :totalcost) for i in idx],
            :counts => [getproperty(out[i], :counts) for i in idx],
            :wcounts => [getproperty(out[i], :wcounts) for i in idx],
            :cweights => [getproperty(out[i], :cweights) for i in idx]
        )
    KMeansResults(assn, costs, other, centers)
end
