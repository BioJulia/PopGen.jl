---
id: tsne
title: t-SNE
sidebar_label: t-SNE
---

t-distributed stochastic neighbor embedding (_a.k.a._ t-SNE) is a dimensionality reduction technique for visualizing high-dimensional data. It does this by giving each datapoint a location in a two or three-dimensional map by minimizing the Kullback–Leibler divergence between the high and low dimensionality probability distributions with respect to the locations of the points in the map. It models each high-dimensional object by a two- or three-dimensional point so similar objects are appear nearer and dissimilar objects appear further apart.

:::caution careful parameterization
Visual clusters can be seriously influenced by the parameters. For example, parameters can be chosen in such a way to identify clusters in data that has none. So, a good understanding of the parameters for t-SNE is necessary. **Although a useful tool, t-SNE is not commonly used in population genetic analysis.** It has been included in this package as a wrapper for [TSNE.jl](https://github.com/lejon/TSne.jl) due to its utility in other disciplines.
:::

## tnse
```julia
tsne(data::PopData, args...; kwargs...)
```
Perform t-SNE (t-Stochastic Neighbor Embedding) on a PopData object, returning a DataFrame. Converts the
PopData object into a matrix of allele frequencies with missing values replaced with
the global mean frequency of that allele. First performs PCA on that matrix, retaining
`reduce_dims` dimensions of the PCA prior to t-SNE analysis. The other positional and keyword arguments
are the same as `tsne` from `TSne.jl`.

### Arguments
- `data`: a `PopData` object
- `ndims`: Dimension of the embedded space (default: `2`)
- `reduce_dims` the number of the first dimensions of X PCA to use for t-SNE, if 0, all available dimension are used (default: `0`)
- `max_iter`: Maximum number of iterations for the optimization (default: `1000`)
- `perplexity`: The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results (default: `30`)
### Keyword Arguments (optional)
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


----
## Acknowledgements
This function is a wrapper for the another package, so really, all the brilliance and effort should be credited to the authors of [TNSE.jl](https://github.com/lejon/TSne.jl).
