function dapc(data::PopData, labels::T = [nothing]; classes::Int64) where T<:AbstractVector
    pca_mtx = permutedims(pca(data))
    #pca_mtx = pca(data)

    if labels != [nothing] 
        labs = PooledArray(labels).refs |> Vector{Int64}
    else
       labs = PooledArray(data.meta.population).refs |> Vector{Int64}
    end
        fit(MulticlassLDA, classes, pca_mtx, labs)
end

dapc2(data::PopData) = multiclass_lda_stats(length(populations(data)), permutedims(pca(data)), Vector{Int64}(PooledArray(data.meta.population).refs))