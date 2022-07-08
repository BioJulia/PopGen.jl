using MultivariateStats: MulticlassLDA

function dapc(data::PopData, labels::T = [nothing]) where T<:AbstractVector#, covestimator_between=SimpleCovariance(), covestimator_within=SimpleCovariance())
    if labels != [nothing] 
        labs = _relabel(labels)
    else
        labs = _relabel(data.sampleinfo.population)
    end
    classes = length(unique(labs))
    _pca = pca(data).proj |> permutedims
    #return labs, _pca
    lda = fit(MulticlassLDA, classes, _pca, labs)
    predict(lda, _pca)
end


function subspacedapc(data::PopData, labels::T = [nothing]; classes::Int64, nda::Int64 = 0, covestimator_between=SimpleCovariance(), covestimator_within=SimpleCovariance()) where T<:AbstractVector
    if labels != [nothing] 
        labs = _relabel(labels)
    else
        labs = _relabel(data.sampleinfo.population)
    end

    lablen = length(unique(labs))
    classes > lablen && error("Number of classes ($classes) must be equal or fewer than number of unique labels ($lablen)")
    
    _pca = pca(data)
    pca_mtx = permutedims(_pca.proj)
    ndakw = iszero(nda) ? min(size(pca_mtx, 1), size(pca_mtx,2) - 1) : nda
    MultivariateStats.fit(SubspaceLDA, classes, pca_mtx, labs, outdim = ndakw, covestimator_between= covestimator_between, covestimator_within = covestimator_within)
end