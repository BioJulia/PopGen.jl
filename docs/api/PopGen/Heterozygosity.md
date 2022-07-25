---
id: heterozygosity
title: Heterozygosity.jl
sidebar_label: Heterozygosity.jl
---
import useBaseUrl from "@docusaurus/useBaseUrl";

<link rel="stylesheet" href={useBaseUrl("katex/katex.min.css")} />

## PopGen.jl/src/Heterozygosity.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 📦 counthet
```julia
counthet(geno::T, allele::Int) where T<:GenoArray
counthet(geno::T, allele::AbstractVector{U}) where T<:GenoArray where U<:Integer
```
Given a `GenoArray`, count the number of times `allele` appears in the
heterozygous state.

----
### 📦 counthom
```julia
counthom(geno::T, allele::Int) where T<:GenoArray
counthom(geno::T, allele::AbstractVector{U}) where T<:GenoArray where U<:Integer
```
Given a `GenoArray`, count the number of times `allele` appears in the
homozygous state.

----
### 📦 _genediversitynei87
```julia
_genediversitynei87(het_exp::T, het_obs::T, n::Union{Integer,T}; corr::Bool = true) where T<: AbstractFloat
_genediversitynei87(het_exp::AbstractFloat, het_obs::Missing, n::Union{Integer,AbstractFloat}; corr::Bool = true)
_genediversitynei87(het_exp::Missing, het_obs::AbstractFloat, n::Union{Integer,AbstractFloat}; corr::Bool = true)
_genediversitynei87(het_exp::Missing, het_obs::Missing, n::Union{Integer,AbstractFloat}; corr::Bool = true)
```
Calculate overall gene diversity with the adjustment/correction given by Nei:

$$1-\sum{\bar{p}^2_i + \frac{H_s}{\tilde{n}\times np} - \frac{H_{obs}}{2\tilde{n}\times np}}$$

- $\tilde{n}$ is the number of genotypes for a locus for a population
- $np$ is the number of genotypes of a locus across all populations, i.e. $\sum{\tilde{n}}$
- $\bar{p}^2_i$ is the observed homozygosity of a locus for that population
- $H_s$ is the within population gene diversity given by:
$$H_s = \frac{\tilde{n}}{\tilde{n}-1} \times 1 - (\bar{p}^2_i - \frac{H_{obs}}{2\tilde{n}})$$

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

Use `corr = false` to ignore sample-size correction `* n/(n-1)`

----
### 📦 _hetero_obs
```julia
_hetero_obs(data::T) where T <: GenoArray
```
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined
as genotypes returning `true` for `_ishet()`. This is numerically feasible because
`true` values are mathematically represented as `1`, whereas `false` are represented
as `0`.

----
### 📦 _hetero_exp
```julia
_hetero_exp(allele_freqs::Vector{T}) where T <: GenoArray
```
Returns the expected heterozygosity of an array of genotypes,
calculated as 1 - sum of the squared allele frequencies.

----

### 🔵 heterozygosity
```julia
heterozygosity(data::PopData; by::Union{Symbol,String} = "locus")
```
Calculate observed and expected heterozygosity in a `PopData` object. For loci,
heterozygosity is calculated in the Nei fashion, such that heterozygosity is
calculated as the average over heterozygosity per locus per population.

**Modes**
- `"locus"` : heterozygosity per locus (default)
- `"sample"` : heterozygosity per individual/sample
- `"population"`: heterozygosity per population

**Example**
```julia
heterozygosity(@nancycats, by = "population" )
```
----

### 📦 _heterozygosity
```julia
_heterozygosity(data::PopData, ::Val{:locus})
_heterozygosity(data::PopData, ::Val{:sample})
_heterozygosity(data::PopData, ::Val{:population})
```

----

### 🔵 samplehet
```julia
samplehet(data::PopData, individual::String)
```
Calculate the observed heterozygosity for an individual in a `PopData` object.
