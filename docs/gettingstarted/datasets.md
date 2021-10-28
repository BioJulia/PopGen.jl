---
id: datasets
title: Provided datasets
sidebar_label: Provided datasets
---

PopGen.jl (via PopGenCore.jl) provides two datasets as examples, `nancycats` and `gulfsharks`. The datasets can be retrieved using the `dataset` function, or their names as macros  (e.g. `@gulfsharks`).

:::info identitcal methods
The methods are identical (one is a wrapper for the other), but the benefit of calling the datasets directly by name is that you get the luxury of tab auto-completion :grin:
:::

## datasets
```julia
PopGen.dataset(::String)
```
Returns a `PopData` object of the dataset you would like to retrieve by calling the dataset as a string by name.

**Example:**
```julia
sharks = PopGen.dataset("gulfsharks")
cats = PopGen.dataset("nancycats")
```
### nancycats

We include the familiar nancycats microsatellite data, as featured in [adegenet](http://adegenet.r-forge.r-project.org), for easy importing into PopGen.jl as `PopData`. As an alternative to `datasets`, you can invoke the `@nancycats` macro.

```julia
julia> cats = @nancycats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
```

The spatial coordinates provided for the dataset in `adegenet` are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted. If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!

### gulfsharks

We also include the SNP dataset used in [Dimens *et al.* 2019](https://link.springer.com/article/10.1007/s00227-019-3533-1) since it was already on hand. Like `nancycats`, we provide a convenient function to load these data into PopGen.jl as `PopData`. As an alternative to `dataset`, you can invoke the `@gulfsharks` macro. 

```julia
julia> sharks = @gulfsharks
PopData{Diploid, 2209 SNP loci}
  Samples: 212
  Populations: 7
  Other Info: ["longitude", "latitude"]
```

