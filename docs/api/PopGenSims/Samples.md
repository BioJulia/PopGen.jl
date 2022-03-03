---
id: popgensims_samples
title: Samples.jl
sidebar_label: Samples.jl
---
## PopGenSims.jl/src/Samples.jl
📦  => not exported | 
🟪 => exported by PopGenSims.jl

### 📦 sample_locus
```julia
sample_locus(locus::Dict, n::Int, ploidy::Signed)
```
Internal function used by `simulate` to take a `Dict` of alleles => frequencies of a locus and return
`n` number of genotypes (n_alleles = `ploidy`) by using weighted sampling of the
allele-frequency pairs.

**Example**
```julia
d = Dict(133 => 0.125,135 => 0.5625,143 => 0.25,137 => 0.0625)

julia> sample_locus(d, 3, 2)
5-element Array{Tuple{Int16,Int16},1}:
 (133, 135)
 (135, 135)
 (143, 137)

julia> sample_locus(d, 3, 3)
5-element Array{Tuple{Int16,Int16,Int16},1}:
 (135, 135, 133)
 (143, 135, 133)
 (137, 135, 135)
```
-----

### 📦 _simulatearbitrary
```julia
_simulatearbitrary(data::PopData, n::Dict{String, Int})
```
-----
### 📦 _simulateflat
```julia
_simulateflat(data::PopData, n::Int)
```
-----

### 📦 _simulatescale
```julia
_simulatescale(data::PopData, scale::Int)
```
-----

### 🟪 simulate
```julia
simulate(data::PopData; n::Int)
```
Simulate `n` number of individuals per population using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object with `n` * `n_populations` samples.

**Example**
```julia
julia> cats = @nanycats;

julia> sims = simulate(cats, n = 100)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 1700
  Populations: 17
```  
----

```julia
simulate(data::PopData; n::Dict{String, Int})
```
Simulate an arbitrary number of samples per populations specified in the Dict `n`, given by `Population => #samples`. Uses
per-population allele frequencies derived from a `PopData` object. 
```julia
julia> cats = @nancycats;

julia> simscheme = Dict("1" => 5, "8" => 3, "11" => 20) ;

julia> simulate(cats, n = simscheme)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 28
  Populations: 3
```
----
```julia
simulate(data::PopData; scale::Int)
```
Simulate individuals per population in the same proportions they appear in the PopData
using per-population allele frequencies. Simulation volume can be multiplied using `scale`,
i.e. if you want to keep the same proportions but generate twice the number of samples, `scale`
would be `2`. Returns a new `PopData` object with `n_samples` * `scale` samples.    

```julia
julia> cats = @nanycats;

julia> sims_prop = simulate(cats, scale = 3)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 711
  Populations: 17
```
