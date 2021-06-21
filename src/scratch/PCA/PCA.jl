using LinearAlgebra, MultivariateStats, CairoMakie, Clustering
include("AlleleMatrices.jl")

# pratio The ratio of variances preserved in the principal subspace.
# mean: data has already been centered
function pca(data::PopData; missings::String = "mean", retain::Float64 = NaN, raw::Bool = false, scale::Bool = false, center::Bool = true)
    mtx = allele_matrix(data, missings = missings, scale = scale, center = center) |> permutedims
    #dims = iszero(dimensions) ? size(mtx, 2) : dimensions
    if isnan(retain)
        M = fit(PCA, mtx, mean = 0, method = :svd)
    else
        pr = retain > 1 ? retain * 0.01 : retain
        M = fit(PCA, mtx, mean = 0, method = :svd, pratio = pr)
    end
    @info "\nRetaining " * string(size(M.proj,2)) * " principal components, capturing ~" * string(round(tprincipalvar(M)/tvar(M) * 100, digits = 1)) * "% of the variance"
    if raw 
        return M
    else
        return MultivariateStats.transform(M, mtx) |> permutedims
    end
end


b = pca(cats)

Makie.scatter(b[:,2], b[:,3], color = colors)