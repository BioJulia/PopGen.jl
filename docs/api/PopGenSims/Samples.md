---
id: popgensims_samples
title: Samples.jl
sidebar_label: Samples.jl
---
PopGenSims.jl/src/Samples.jl
❗ => not exported | 
⚫ => exported by PopGenSims.jl

### ❗sample_locus
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

### ⚫ simulate
    simulate(data::PopData; n::Int = 100)
Simulate `n` number of individuals (default: `100`) per population using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object.

**Example**
```julia
cats = @nancycats;

julia> sims = simulate(cats , n = 100)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 1700
  Populations: 17
  
julia> sampleinfo(sims)

  1700×5 DataFrame
  Row │ name      population  ploidy   
      │ String    String      Int8      
──────┼───────────────────────────────
    1 │ sim_1     1                2    
    2 │ sim_2     1                2    
    3 │ sim_3     1                2    
    4 │ sim_4     1                2    
    5 │ sim_5     1                2    
  ⋮   │    ⋮          ⋮         ⋮ 
 1697 │ sim_1697  17               2  
 1698 │ sim_1698  17               2  
 1699 │ sim_1699  17               2  
 1700 │ sim_1700  17               2  
                                         1691 rows omitted

julia> genodata(sims)
15300×4 DataFrame
   Row │ name      population  locus   genotype   
       │ String    String      String  Tuple…?    
───────┼──────────────────────────────────────────
     1 │ sim_1     1           fca8    (135, 143)
     2 │ sim_1     1           fca23   (136, 146)
     3 │ sim_1     1           fca43   (141, 145)
     4 │ sim_1     1           fca45   (120, 126)
     5 │ sim_1     1           fca77   (156, 156)
   ⋮   │    ⋮          ⋮         ⋮         ⋮
 15297 │ sim_1700  17          fca78   (150, 150)
 15298 │ sim_1700  17          fca90   (197, 197)
 15299 │ sim_1700  17          fca96   (113, 113)
 15300 │ sim_1700  17          fca37   (208, 208)
                                15291 rows omitted
```
```