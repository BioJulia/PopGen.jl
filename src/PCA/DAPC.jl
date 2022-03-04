function pca(data::PopData; maxpc::Int = 0, method::Symbol = :svd, pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
    meankw = center ? 0 : nothing
    mtx = _allelematrix(data, scale = scale, center = center)
    pckw = iszero(maxpc) ? min(size(mtx, 1), size(mtx,2) - 1) : maxpc
    MultivariateStats.fit(PCA, mtx; maxoutdim=pckw, mean = meankw, pratio = pratio, method = method)
end


function dapc(data::PopData, labels::T = [nothing]; classes::Int64, nda::Int64 = 0, covestimator_between=SimpleCovariance(), covestimator_within=SimpleCovariance()) where T<:AbstractVector
    if labels != [nothing] 
        labs = _relabel(labels)
    else
        labs = _relabel(data.sampleinfo.population)
    end

    lablen = length(unique(labs))
    classes > lablen && error("Number of classes ($classes) must be equal or fewer than number of unique labels ($lablen)")
    
    _pca = pca(data)
    return _pca
    pca_mtx = permutedims(_pca.proj)
    ndakw = iszero(nda) ? min(size(pca_mtx, 1), size(pca_mtx,2) - 1) : nda
    MultivariateStats.fit(MulticlassLDA, classes, pca_mtx, labs, outdim = ndakw, covestimator_between= covestimator_between, covestimator_within = covestimator_within)
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