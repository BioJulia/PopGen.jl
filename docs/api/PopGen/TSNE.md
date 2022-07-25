---
id: tsne
title: TSNE.jl
sidebar_label: TSNE.jl
---
## PopGen.jl/src/SummaryInfo.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 🔵 tsne
```julia
tsne(data::PopData, args...; kwargs...)
```
Perform t-SNE (t-Stochastic Neighbor Embedding) on a PopData object, returning a DataFrame. Converts the
PopData object into a matrix of allele frequencies with missing values replaced with
the global mean frequency of that allele. First performs PCA on that matrix, retaining
`reduce_dims` dimensions of the PCA prior to t-SNE analysis. The positional and keyword arguments
are the same as `tsne` from `TSne.jl`.

**Arguments**
- `data`: a `PopData` object
- `ndims`: Dimension of the embedded space (default: `2`)
- `reduce_dims` the number of the first dimensions of X PCA to use for t-SNE, if 0, all available dimension are used (default: `0`)
- `max_iter`: Maximum number of iterations for the optimization (default: `1000`)
- `perplexity`: The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results (default: `30`)

**Keyword Arguments (optional)**
- `distance`: type `Function` or `Distances.SemiMetric`, specifies the function to
  use for calculating the distances between the rows
- `pca_init`: whether to use the first `ndims` of the PCA as the initial t-SNE layout,
  if `false` (the default), the method is initialized with the random layout
- `max_iter`: how many iterations of t-SNE to do
- `verbose`: output informational and diagnostic messages
- `progress`: display progress meter during t-SNE optimization (default: `true`)
- `min_gain`: `eta`: `initial_momentum`, `final_momentum`, `momentum_switch_iter`,
  `stop_cheat_iter`: `cheat_scale` low-level parameters of t-SNE optimization
- `extended_output`: if `true`, returns a tuple of embedded coordinates matrix,
  point perplexities and final Kullback-Leibler divergence