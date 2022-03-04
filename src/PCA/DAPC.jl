function pca(data::PopData; maxpc::Int = 0, method::Symbol = :svd, pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
    meankw = center ? 0 : nothing
    mtx = _allelematrix(data, scale = scale, center = center)
    pckw = iszero(maxpc) ? min(size(mtx, 1), size(mtx,2) - 1) : maxpc
    MultivariateStats.fit(PCA, mtx; maxoutdim=pckw, mean = meankw, pratio = pratio, method = method)
end


function dapc(data::PopData, labels::T = [nothing]; classes::Int64) where T<:AbstractVector
    if labels != [nothing] 
        labs = PooledArray(labels).refs |> Vector{Int64}
    else
       labs = PooledArray(data.sampleinfo.population).refs |> Vector{Int64}
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