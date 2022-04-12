"""
    pca(data::PopData; maxpc::Int = 0, method::Symbol = :svd, missings::String = "mean", pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)

Perform a Principal Component Analysis on a PopData object. Returns an indexible `MultivariateStats.PCA` object.

#### Arguments
- `data::PopData`: a `PopData` object
#### keyword arguments
- `method::Symbol`: The PCA method to use (default: `:svd`)
    - `:cov`: based on covariance matrix decomposition
    - `:svd`: based on Singular Value Decomposition of the input data
- `maxpc::Int`: The maximum number of principal components to retain (default: 0 = `(min(d, ncol-1))`)
- `missings::String`: How to treat missing genotypes in the allele frequency matrix (default: `mean`)
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
    - `"missing"`: keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
- `pratio::Float64`: The maxium ratio of variances preserved in the principal subspace (default = `0.99`)
- `center::Bool`: whether to center the allele frequency matrix (default: `false`)
- `scale::Bool`: whether to Z-score scale the allele frequency matrix (default: `true`)
"""
function pca(data::PopData; maxpc::Int = 0, method::Symbol = :svd, missings::String = "mean", pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
    meankw = center ? 0 : nothing
    mtx = _allelematrix(data, scale = scale, center = center, missings = missings)
    pckw = iszero(maxpc) ? min(size(mtx, 1), size(mtx,2) - 1) : maxpc
    fit(PCA, mtx; maxoutdim=pckw, mean = meankw, pratio = pratio, method = method)
end
