---
id: pca
title: Principal Component Analysis
sidebar_label: Principal Component Analysis
---

A common way to analyze genetic data is dimensionality reduction, and PopGen.jl provides `pca()` to perform a Pricipal Component Analysis. This
function wraps `fit(PCA, ...)` from ([MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl)) to be used on `PopData` objects.
The genotypes are processed into a matrix of (rows: samples, cols: allele frequencies), giving you the option of how to manage `missing` data, and the PCA is performed on this allele-frequency matrix. 

:::tip suppressing output
For datasets greater than 10 loci, we recommend appending a semicolon to the end of the function call to suppress output to the REPL. ([issue #186](https://github.com/JuliaStats/MultivariateStats.jl/issues/186))
:::

### Principal Component Anaylsis
```julia
pca(::PopData; maxpc::Int = 0, method::Symbol = :svd, missings::String = "mean", pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
```
Perform a Principal Component Analysis on a PopData object. Returns an indexible `MultivariateStats.PCA` object.

#### keyword arguments
- `method::Symbol`: The PCA method to use (default: `:svd`)
    - `:cov`: based on covariance matrix decomposition
    - `:svd`: based on Singular Value Decomposition of the input data
- `maxpc::Int`: The maximum number of principal components to retain (default: 0 = `(min(d, ncol-1))`)
- `missings::String`: How to treat missing genotypes in the allele frequency matrix (default: `"mean"`)
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
    - `"missing"`: keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
- `pratio::Float64`: The maxium ratio of variances preserved in the principal subspace (default = `0.99`)
- `center::Bool`: whether to center the allele frequency matrix (default: `false`)
- `scale::Bool`: whether to Z-score scale the allele frequency matrix (default: `true`)

**Example**
```julia
julia> cats = @nancycats;

julia> pca_cats = pca(cats, maxpc = 40); 

julia> pca_cats.proj       # the projection matrix
237×40 Matrix{Float64}:
 -0.0922663  -0.00353549  -0.0205698    …   0.0148947    0.0114937     0.0119581    
 -0.0939368   0.00935327  -0.0488424       -0.0167817   -0.0144753    -0.0348131    
 -0.0545479  -0.00597291   0.0527322       -0.0359343    0.00588278    0.0205089    
 -0.0745392   0.0306307   -0.0210157        0.100871     0.0440558     0.0728491    
  ⋮                                     ⋱
 -0.0859615   0.0369075    0.000597856     -0.0633689    0.0249073    -0.0405845    
 -0.0547981   0.0322698   -0.0107031       -0.1138      -0.000545308   0.00296074   
 -0.0950973   0.0163225    0.00588324   …   0.0256939    0.0831163     0.00374728   
 -0.0843321  -0.0082427   -0.0309442        0.00382228   0.0179109     0.0217293

julia> pca_cats.prinvars   # the principal variances
40-element Vector{Float64}:
 20.906886724195086
  8.015966142401036
  7.501281679600134
  6.628421860444929
  ⋮
  2.6556826537989755
  2.599100605083436
  2.5912046561168447
  2.5000331618590637
```

---------------------
## Acknowledgements
The PCA method is an extension of the hyper-fast methods written in [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl). The centering and scaling are also outsourced to the methods present in [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl). A lot of clever package authors and contributors in this ecosystem!