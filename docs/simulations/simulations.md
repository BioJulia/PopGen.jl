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
simulate(data::PopData; n::Int)
simulate(data::PopData; n::Dict{String,Int})
simulate(data::PopData; scale::Int)
```
Simulate data using per-population allele frequencies. The simulated samples will have the naming convention `sim_1...sim_#` where `#` is the total number of simulated samples in the new PopData.

<Tabs
  block={true}
  defaultValue="f"
  values={[
    { label: 'fixed samples', value: 'f', },
    { label: 'arbitrary samples', value: 'a', },
    { label: 'proportional samples', value: 'p', },
  ]
}>
<TabItem value="f">

```julia
simulate(data::PopData; n::Int)
```
If you want to simulate data with a fixed number of individuals per population, you can do so with `simulate(PopData, n = Int)`, which takes a `PopData` object and simulates `n` number of individuals per population. Returns a new PopData with samples having the same ploidy as the source `PopData`, but will **not** work on mixed-ploidy data. 

In the example below, we simulate 100 individuals per
population using the nancycats data, which has 17 populations, therefore the resulting `PopData` will have 1700 samples (100 samples x 17 populations).

**Example**
```julia
julia> cats = @nancycats;

julia> sims = simulate(cats , n = 100)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 1700
  Populations: 17
```

</TabItem>
<TabItem value="a">

```julia
simulate(data::PopData; n::Dict{Population,Int})
```
If you want to simulate an arbitrary number of individuals for arbitrary populations, use `simulate(PopData, n = Dict{String, Int})`, which takes a `PopData` object and simulates samples within populations as specified in the input `Dict`. Returns a new PopData with samples having the same ploidy as the source `PopData`, but will **not** work on mixed-ploidy data.

In the example below, we create a dictionary with the notation `Population => #samples` to simulate a specific number of samples for 3 particular populations. The resulting PopData will have 28 samples (5+3+20) across 3 populations ("1", "8", "11").

```julia
julia> cats = @nancycats;

julia> simscheme = Dict("1" => 5, "8" => 3, "11" => 20) ;

julia> simulate(cats, n = simscheme)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 28
  Populations: 3
```

</TabItem>
<TabItem value="p">

```julia
simulate(data::PopData; n::Int)
```
If you want to simulate data while keeping the proportion of individuals per population consistent with the source PopData, use `simulate(PopData, scale = Int)`, which takes a `PopData` object and simulates the same number of individuals per population multiplied by `scale` (i.e. if `scale=2`, there will be twice the number of simulated individuals compared to the original PopData). Returns a new PopData with samples having the same ploidy as the source `PopData`, but will **not** work on mixed-ploidy data. 

In the example below, we simulate 3x the number of samples of the original nancycats data, which has 237 samples x 17 populations, therefore the resulting `PopData` will have 711 samples (237 samples x 3). In this example, each population will have 3x the number of samples as the original nancycats data.

**Example**
```julia
julia> cats = @nancycats;

julia> sims = simulate(cats , scale = 1)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 711
  Populations: 17
```

</TabItem>
</Tabs>