include("AlleleMatrices.jl")

# pratio The ratio of variances preserved in the principal subspace.
# mean: data has already been centered
function pca(data::PopData; missings::String = "mean", dimensions::Int = 0, var_ratio::Float64 = 0.0)
    mtx = allele_matrix(data, missings = "mean", scale = false, center = true) |> permutedims
    dims = iszero(dimensions) ? size(mtx, 2) : dimensions
    if iszero(var_ratio)
        M = fit(PCA, mtx, maxoutdim=dims, mean = 0)
    else
        M = fit(PCA, mtx, maxoutdim=dims, mean = 0, pratio = var_ratio)
    end
    @info "\nSelected " * string(size(M.proj,2)) * " dimensions capturing " * string(tprincipalvar(M) / tvar(M) * 100) * "% of the variance"
    MultivariateStats.transform(M, mtx) |> permutedims
end

a = pca(sharks)
b = pca(sharks)

p2 = scatter(b[:,2], b[:,3], xlabel="PC1", c = colors, ylabel="PC2", title="PCA")

# Columns correspond to observations

X = allele_matrix(cats, missings = "mean", scale = true, center = true)

# Fit PCA model in 2d, for visualization
M = fit(PCA, X'; maxoutdim=2)
             
# Transform data
x = MultivariateStats.transform(M, X')

# Plot dims 1 and 2 w/o PCA
p1 = scatter(X[1,:], X[2,:], xlabel="Dim. 1", ylabel="Dim. 2", title="Original")

# Plot dims 1 and 2 with PCA
p2 = scatter(x[1,:], x[2,:], xlabel="PC1", c = colors, ylabel="PC2", title="PCA")

colors = deepcopy(cats.meta.population)
replace!(colors, 
    "1" => "antiquewhite3",
    "2" => "antiquewhite4",
    "3" => "aqua"         ,
    "4" => "aquamarine"   ,
    "5" => "aquamarine1"  ,
    "6" => "aquamarine2"  ,
    "7" => "aquamarine3"  ,
    "8" => "aquamarine4"  ,
    "9" => "azure"        ,
    "10" => "azure1"      , 
    "11" => "azure2"      , 
    "12" => "azure3"      , 
    "13" => "azure4"      , 
    "14" => "beige"       , 
    "15" => "bisque"      , 
    "16" => "bisque1"     , 
    "17" => "bisque2"
)


colors = deepcopy(sharks.meta.population)
replace!(colors,
    "Cape Canaveral" => "#0a2463", 
    "Georgia"        => "#3e92cc", 
    "South Carolina" => "#fffaff", 
    "Florida Keys"   => "#d8315b", 
    "Mideast Gulf"   => "#1e1b18",
    "Northeast Gulf" => "#78ffd6",
    "Southeast Gulf" => "#b33f62"
)

# Compare plots
plot(p1, p2, layout=2, legend=false)

# manual pca
x = X
#x = x .- mean(x, dims=1)     # Centre the data
xcov = Symmetric(cov(x))    #Compute the covariance matrix
ef = eigen(xcov)    #Get eigenvalues and eigenvectors
eigvalvec = reverse(ef.values, dims=1)    #Orient eigenvalues largest to smallest
eigvecmat = reverse(ef.vectors, dims=2)    #Orient eigenvectors the same
eigvalcs = cumsum(eigvalvec)    #Cumulative sum of eigenvalues
vratio = eigvalcs ./ eigvalcs[end]    #Ratios of variance explained
pc = x * eigvecmat     #Principal components ordered largest eigenvalue to smallest
#
p2 = scatter(pc[:,1], pc[:,3], xlabel="PC1", c = colors, ylabel="PC2", title="PCA")
