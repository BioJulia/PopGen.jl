function kmeans(data::PopData; k_range::UnitRange{Int64} = 0:0, iterations::Int64 = 100)
    kr = k_range == 0:0 ? (2:(size(data)[1] ÷ 10)) : k_range
    pc = pca(data) |> permutedims
    map(i -> Clustering.kmeans(pc, i, maxiter = iterations), kr)
end