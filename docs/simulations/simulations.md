---
id: simulate_samples
title: Simulating Samples
sidebar_label: Simulating Samples
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

:::note Requires PopGenSims.jl
To perfom simulations, you will need add and import the package `PopGenSims.jl` (available [here](https://github.com/pdimens/PopGenSims.jl)).
:::

## Simulate samples within populations
```julia
simulate(data::PopData; n::Int = 100)
```
If you want to generate simulated data of a certain number of individuals per population, you can do so with the `simulate()` function, which takes a `PopData` object and simulates `n` number of individuals per population using the allele frequencies of each population. This returns 
new `PopData`. The simulated samples will have the naming convention `sim_#` where `#` is a number from 1:`n`. These simulations return samples with the same ploidy as the source `PopData`, but will **not** work on mixed-ploidy data. 

In the example below, we simulate 100 individuals per
population using the nancycats data, which has 17 populations, therefore the resulting `PopData` will have 1700 samples (100 samples x 17 populations)

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
```

Here is a look inside the `PopData` to verify everything looks as expected.

<Tabs
  block={true}
  defaultValue="m"
  values={[
    { label: 'meta', value: 'm', },
    { label: 'loci', value: 'l', },
  ]
}>
<TabItem value="m">

```
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
```

</TabItem>
<TabItem value="l">

```
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

</TabItem>
</Tabs>