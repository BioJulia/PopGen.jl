---
id: tsne
title: t-SNE
sidebar_label: t-SNE
---

:::note
PopGen.jl v0.10x removed TSNE.jl as a dependency and instead this page is inteded to walk you through how to use PopGen.jl and [TSne.jl](https://github.com/lejon/TSne.jl) together. You will need to install TSNE.jl with
```julia
julia> ]add TSNE
# or #
julia> using Pkg; Pkg.add("TSne")
```
:::

t-distributed stochastic neighbor embedding (_a.k.a._ t-SNE) is a dimensionality reduction technique for visualizing high-dimensional data. It does this by giving each datapoint a location in a two or three-dimensional map by minimizing the Kullback–Leibler divergence between the high and low dimensionality probability distributions with respect to the locations of the points in the map. It models each high-dimensional object by a two- or three-dimensional point so similar objects are appear nearer and dissimilar objects appear further apart. Read more about it [here](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding).

:::caution careful parameterization
Visual clusters can be seriously influenced by the parameters. For example, parameters can be chosen in such a way to identify clusters in data that has none. So, a good understanding of the parameters for t-SNE is necessary. **Although a useful tool, t-SNE is not commonly used in population genetic analysis.** It is included here due to its utility in adjacent disciplines.
:::

## Step 1: Create the allele frequency matrix
You will need to input an allele frequency matrix into the `tsne` function. PopGenCore.jl has a simple function to do that, but it's not exported. We will assume your `PopData` object is called `data` in this example. Missing values in the matrix will be replaced with the global mean frequency of that allele.
```julia
mtx = PopGenCore.freqmatrix_mean(data)
```

## Step 2: Perform t-SNE on the frequency matrix
Once you have the allele frequency matrix, you just plug it into `tsne` from `TSne.jl` as the first positional argument.
```julia
results_tsne = tsne(mtx, ndim, reduce_dims, max_iter, perplexity; kwargs...)
```
### Arguments
- `ndim`: Dimension of the embedded space (default: `2`)
- `reduce_dims` the number of the first dimensions of X PCA to use for t-SNE, if 0, all available dimension are used (default: `0`)
- `max_iter`: how many iterations of t-SNE to do (default: `1000`)
- `perplexity`: the perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results (default: `30`)

### Keyword Arguments (optional)
- `distance`: if `true`, specifies that `X` is a distance matrix, if of type `Function` or `Distances.SemiMetric`, specifies the function to use for calculating the distances between the rows (or elements, if `X` is a vector) of `X`
- `pca_init`: whether to use the first ndims of X PCA as the initial t-SNE layout, if `false` (the default), the method is initialized with the random layout
- `verbose` if `true`, output informational and diagnostic messages
- `progress`: if `true`, display progress meter during t-SNE optimization
- `min_gain`, `eta`, `initial_momentum`, `final_momentum`, `momentum_switch_iter`, `stop_cheat_iter`, `cheat_scale`: low-level parameters of t-SNE optimization
- `extended_output`: if `true`, returns a tuple of embedded coordinates matrix, point perplexities and final Kullback-Leibler divergence



## Step 3: Nicer formatting for results
You may want to make the results look nicer and more legible, which can be done with this code snippet. You may need to import `DataFrames`, since PopGen.jl doesn't reexport the `DataFrame` function.
```julia
df_tsne = DataFrame(results_tsne, String["dim"*"$i" for i in 1:size(results_tsne, 2)])
insertcols!(df_tsne, 1, :name => data.sampleinfo.name, :population => data.sampleinfo.population)
```

