sharks = @gulfsharks

function dapc(data::PopData, labels::T = [nothing]; classes::Int64) where T<:AbstractVector
    if labels != [nothing] 
        labs = PooledArray(labels).refs |> Vector{Int64}
    else
       labs = PooledArray(data.meta.population).refs |> Vector{Int64}
    end

    lablen = length(unique(labs))
    classes < lablen && error("Number of classes ($classes) must be equal or fewer than number of unique labels ($lablen)")
    
    pca_mtx = permutedims(pca(data))
    #pca_mtx = pca(data)
    #pca_mtx = allele_matrix(data, missings = "mean", scale = true, center = false)

    fit(MulticlassLDA, classes, pca_mtx, labs)
end

tst = dapc(sharks, classes = 9)

dapc2(data::PopData) = multiclass_lda_stats(length(populations(data)), permutedims(pca(data)), Vector{Int64}(PooledArray(data.meta.population).refs))