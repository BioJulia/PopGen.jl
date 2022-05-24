---
id: pairwisekinship
title: PairwiseKinship.jl
sidebar_label: PairwiseKinship.jl
---
## PopGen.jl/src/Kinship/PairwiseKinship.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 📦 _bootstrapsummary
```julia
_bootstrapsummary(::Vector{Union{Missing, Float64}}, width::Tuple{Float64, Float64})
```
Return the mean, median, standard error, and quantiles (given by `witdth`) of relatedness resampling.


----

### 📦 _bootstrapgenos_all
```julia
_bootstrapgenos_all(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::NamedTuple; method::Function, iterations::Int)
```
Perform `iterations` number of bootstrap resampling iterations of all genotypes between pair (`ind1` `ind2`). Returns a vector of length `interatotions`
of the relatedness estimate given by method `method`. This is an internal function with `locus_names`, `n_per_loc`, and `alleles` supplied by `relatedness_boot_all`.

----

### 📦 _bootstrapgenos_nonmissing
```julia
bootstrapgenos_nonmissing(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}, n_per_loc::Vector{Int}, alleles::NamedTuple; method::Function, iterations::Int)
```
Perform `iterations` number of bootstrap resampling iterations of only shared (nonmissing) genotypes between pair (`ind1` `ind2`). Returns a vector of length `interatotions`
of the relatedness estimate given by method `method`. This is an internal function with `locus_names`, `n_per_loc`, and `alleles` supplied by `relatedness_boot_nonmissing`.


----

### 📦 _kinship_boot_all
```julia
_kinship_boot_all(::PopData, sample_names::Vector{String}; method::Function, iterations::Int, interval::Tuple{Float64, Float64})
```
Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. Bootstrapping resamples using
the `all` method, where resampling occurs over all loci. This is an internal function with all arguments provided by `relatedness`.


----

### 📦 _kinship_boot_nonmissing
```julia
_kinship_boot_nonmissing(::PopData, sample_names::Vector{String}; method::F, iterations::Int, interval::Tuple{Float64, Float64}) where F
```
Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. Bootstrapping resamples using
the `nonmissing` method, where resampling occurs over only shared non-missing loci. This is an internal function with all arguments provided by `relatedness`.


----

### 📦 _kinship_noboot
```julia
_kinship_noboot(::PopData, sample_names::Vector{String}; method::F) where F
```
Calculate pairwise relatedness between all combinations of the provided `sample_names` for each `method` provided. 
This is an internal function with arguments provided by `relatedness`.


----

### 🔵 kinship
```julia
# compare all samples
kinship(::PopData; method::Function, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String, inbreeding::Bool = false)
# to compare specific samples
kinship(::PopData, samples; method::F, iterations::Int64, interval::Tuple{Float64, Float64}, resample::String, inbreeding::Bool = false)
```
Return a dataframe of pairwise relatedness estimates for all or select pairs of `samples` in a `PopData` object using 
method(s) `F` where `F` is one or several of the methods listed below. If no bootstrapping is required, then the only 
necessary keyword to provide is `method = ` and `inbreeding = ` for the `dyadicLikelihood` method (see examples below). **Note:** samples must be diploid.

#### Estimator methods

The available estimators are listed below and are functions themselves. `relatedness` takes the
function names as arguments (**case sensitive**), therefore do not use quotes or colons
in specifying the methods. Multiple methods can be supplied as a vector. All of these methods will tab-autocomplete.
For more information on a specific method, please see the respective docstring (e.g. `?Loiselle`).
- `Blouin`
- `dyadicLikelihood`
- `LiHorvitz`
- `Loiselle`
- `Lynch`
- `LynchLi`
- `LynchRitland`
- `Moran`
- `QuellerGoodnight`
- `Ritland`
- `Wang`

#### Inbreeding

Use the `inbreeding` keyword to specify whether to allow inbreeding (`true`) or not (`false`, default).
This is only relevant for the `dyadicLikelihood` method.

#### Bootstrapping
To calculate means, medians, standard errors, and confidence intervals using bootstrapping,
set `iterations = n` where `n` is an integer greater than `0` (the default) corresponding to the number
of bootstrap iterations you wish to perform for each pair. The default confidence interval is `(0.05, 0.95)` (i.e. 90%),
however that can be changed by supplying the keyword `interval = (low, high)` where `low` and `high` are the intervals you want 
(as `AbstractFloat`). The returned DataFrame will have 5 columns per `method` with bootstrapped parameters having the naming
convention of `Method_parameter`. The output may have more columns than will fit on your screen, so `DataFrames.names(out_df)`
may be useful to see a list of the column names.

#### Resampling methods

There are two available resampling methods, `"all"` (default  & recommended) and `"nonmissing"`.
- `"all"` : resamples all loci for a pair of individuals and then drops missing loci between them
    - speed: slower
    - pro: better resampling variation
    - con: by chance some iterations may have a lot of missing loci that have to be dropped
- `"nonmissing"` : resamples only the shared non-missing loci between the pair
    - speed: faster
    - pro: every iteration guarantees the same number of loci compared between the pair
    - con: too-tight confidence intervals due to less possible variation

#### Examples
```
julia> cats = @nancycats;

julia> kinship(cats, method = Ritland)
27966×4 DataFrame
│ Row   │ sample_1 │ sample_2 │ n_loci │ Ritland    │
│       │ String   │ String   │ Int64  │ Float64?   │
├───────┼──────────┼──────────┼────────┼────────────┤
│ 1     │ N215     │ N216     │ 8      │ 0.258824   │
│ 2     │ N215     │ N217     │ 8      │ 0.193238   │
│ 3     │ N215     │ N218     │ 8      │ 0.127497   │
⋮
│ 27964 │ N281     │ N289     │ 8      │ 0.0892068  │
│ 27965 │ N281     │ N290     │ 7      │ 0.104614   │
│ 27966 │ N289     │ N290     │ 7      │ 0.0511663  │

julia> kinship(cats, ["N7", "N111", "N115"], method = [Ritland, Wang])
3×5 DataFrame
│ Row │ sample_1 │ sample_2 │ n_loci │ Ritland    │ Wang      │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?  │
├─────┼──────────┼──────────┼────────┼────────────┼───────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.129432  │ -0.395806 │
│ 2   │ N7       │ N115     │ 9      │ -0.0183925 │ 0.0024775 │
│ 3   │ N111     │ N115     │ 9      │ 0.0240152  │ 0.183966  │

julia> kinship(cats, ["N7", "N111", "N115"], method = [Loiselle, Moran], iterations = 100, interval = (0.025, 0.975))
3×13 DataFrame. Omitted printing of 7 columns
│ Row │ sample_1 │ sample_2 │ n_loci │ Loiselle   │ Loiselle_mean │ Loiselle_median │
│     │ String   │ String   │ Int64  │ Float64?   │ Float64?      │ Float64?        │
├─────┼──────────┼──────────┼────────┼────────────┼───────────────┼─────────────────┤
│ 1   │ N7       │ N111     │ 9      │ -0.101618  │ 0.0141364     │ 0.00703954      │
│ 2   │ N7       │ N115     │ 9      │ -0.0428898 │ 0.0743497     │ 0.0798708       │
│ 3   │ N111     │ N115     │ 9      │ 0.13681    │ 0.266043      │ 0.257748        │

julia> DataFrames.names(ans)
13-element Array{String,1}:
 "sample_1"
 "sample_2"
 "n_loci"
 "Loiselle"
 "Loiselle_mean"
 "Loiselle_median"
 "Loiselle_SE"
 "Loiselle_CI_95"
 "Moran"
 "Moran_mean"
 "Moran_median"
 "Moran_SE"
 "Moran_CI_95"
```

### 🔵 merge_kinship
```julia
merge_kinship(data::NamedTuple)
```
A convenience function that takes the `NamedTuple` output from `kinship` performed with bootstrapping
and returns one large DataFrame.
