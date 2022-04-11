function pca(data::PopData; maxpc::Int = 0, method::Symbol = :svd, pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
    meankw = center ? 0 : nothing
    mtx = _allelematrix(data, scale = scale, center = center)
    pckw = iszero(maxpc) ? min(size(mtx, 1), size(mtx,2) - 1) : maxpc
    fit(PCA, mtx; maxoutdim=pckw, mean = meankw, pratio = pratio, method = method)
end
