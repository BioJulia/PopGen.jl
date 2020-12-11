---
id: popgensims_samples
title: Samples.jl
sidebar_label: Samples.jl
---

## `simulate`
    simulate(data::PopData; n::Int = 100)
Simulate `n` number of individuals (default: `100`) per population using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object.
#### Example
```julia
cats = @nancycats;

julia> sims = simulate(cats , n = 100)
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 1700
  Loci: 9
  Populations: 17
  Coordinates: absent

julia> sims.meta
  1700×5 DataFrame
  │ Row  │ name     │ population │ ploidy │ longitude │ latitude │
  │      │ String   │ String     │ Int64  │ Missing   │ Missing  │
  ├──────┼──────────┼────────────┼────────┼───────────┼──────────┤
  │ 1    │ sim_1    │ 1          │ 2      │ missing   │ missing  │
  │ 2    │ sim_2    │ 1          │ 2      │ missing   │ missing  │
  │ 3    │ sim_3    │ 1          │ 2      │ missing   │ missing  │
  │ 4    │ sim_4    │ 1          │ 2      │ missing   │ missing  │
  ⋮
  │ 1696 │ sim_1696 │ 17         │ 2      │ missing   │ missing  │
  │ 1697 │ sim_1697 │ 17         │ 2      │ missing   │ missing  │
  │ 1698 │ sim_1698 │ 17         │ 2      │ missing   │ missing  │
  │ 1699 │ sim_1699 │ 17         │ 2      │ missing   │ missing  │
  │ 1700 │ sim_1700 │ 17         │ 2      │ missing   │ missing  │  
  
julia> sims.loci
  15300×4 DataFrame
  │ Row   │ name     │ population │ locus  │ genotype   │
  │       │ String   │ String     │ String │ Tuple…?    │
  ├───────┼──────────┼────────────┼────────┼────────────┤
  │ 1     │ sim_1    │ 1          │ fca8   │ (135, 135) │
  │ 2     │ sim_1    │ 1          │ fca23  │ (132, 140) │
  │ 3     │ sim_1    │ 1          │ fca43  │ (139, 139) │
  │ 4     │ sim_1    │ 1          │ fca45  │ (126, 126) │
  ⋮
  │ 15297 │ sim_1700 │ 17         │ fca78  │ (142, 142) │
  │ 15298 │ sim_1700 │ 17         │ fca90  │ (199, 199) │
  │ 15299 │ sim_1700 │ 17         │ fca96  │ (113, 113) │
  │ 15300 │ sim_1700 │ 17         │ fca37  │ (208, 208) │
```