---
id: heterozygosity
title: Heterozygosity.jl
sidebar_label: Heterozygosity.jl
---
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />

### `ishom`
```julia
ishom(locus::T) where T <: GenotypeArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if it is, `false` if it isn't, and `missing` if it's `missing`. The vector version simply broadcasts the function over the elements.

----

### `ishet`
```julia
ishet(locus::T) where T <: GenotypeArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if it is, `false` if it isn't. The vector version simply broadcasts the function over the elements. Under the hood, this function is simply `!ishom`.

----

### `gene_diversity_nei87`
```julia
gene_diversity_nei87(het_exp::Union{Missing,AbstractFloat}, het_obs::Union{Missing,AbstractFloat}, n::Union{Integer, Float64}, corr::Bool = true)
```
Calculate overall gene diversity with the adjustment/correction given by Nei:

$$1-\sum{\bar{p}^2_i + \frac{H_s}{\tilde{n}\times np} - \frac{H_{obs}}{2\tilde{n}\times np}}$$

- $\tilde{n}$ is the number of genotypes for a locus for a population
- $np$ is the number of genotypes of a locus across all populations, i.e. $\sum{\tilde{n}}$
- $\bar{p}^2_i$ is the observed homozygosity of a locus for that population
- $H_s$ is the within population gene diversity given by:
$$H_s = \frac{\tilde{n}}{\tilde{n}-1} \times 1 - (\bar{p}^2_i - \frac{H_{obs}}{2\tilde{n}})$$

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

Use `corr = false` to ignore sample-size correction `* n/(n-1)`.

----

### `hetero_o`
```julia
hetero_o(data::T) where T <: GenotypeArray
```
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined as genotypes returning `true` for `ishet()`. This is numerically feasible because `true` values are mathematically represented as `1`, whereas `false` are represented as `0`.

----

### `hetero_e`
```julia
hetero_e(allele_freqs::Vector{T}) where T <: GenotypeArray
```
Returns the expected heterozygosity of an array of genotypes, calculated as 1 - sum of the squared allele frequencies.

----

### `heterozygosity`
```julia
heterozygosity(data::PopData, by::String = "locus")
```
Calculate observed and expected heterozygosity in a `PopData` object. For loci, heterozygosity is calculated in the Nei fashion, such that heterozygosity is calculated as the average over heterozygosity per locus per population.

**Modes**
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population

**Example**
```julia
heterozygosity(@nancycats, "population" )
```

----

### `het_sample`
```julia
het_sample(data::PopData, individual::String)
```
Calculate the observed heterozygosity for an individual in a `PopData` object. Returns an array of heterozygosity values.
