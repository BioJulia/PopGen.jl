---
id: pca
title: Principal Component Analysis
sidebar_label: Principal Component Analysis
---

A common way to analyze genetic data is dimensionality reduction, and PopGen.jl provides `pca()` to perform a Pricipal Component Analysis. This
function wraps `fit(PCA, ...)` from `MultivariateStats.jl` ([link](https://github.com/JuliaStats/MultivariateStats.jl)) to be used on `PopData` objects.
The genotypes are processed into a matrix of (rows: samples, cols: allele frequencies), giving you the option of how to manage `missing` data, and the PCA is performed on this allele-frequency matrix. 

:::tip suppressing output
For datasets greater than 10 loci, we recommend appending a semicolon to the end of the function call to suppress output to the REPL. ([issue #186](https://github.com/JuliaStats/MultivariateStats.jl/issues/186))
:::

### `pca`

```julia
pca(::PopData; maxpc::Int = 0, method::Symbol = :svd, missings::String = "mean", pratio::Float64 = 0.99, center::Bool = false, scale::Bool = true)
```
Perform a Principal Component Analysis on a PopData object. Returns an indexible `MultivariateStats.PCA` object.

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

**Example**
```julia
julia> cats = @nancycats;

julia> pca_cats = pca(cats, maxpc = 40); 

julia> pca_cats.proj       # the projection matrix



julia> pca_cats.prinvars  # the principal variances
```