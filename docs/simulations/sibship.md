---
id: sibship_simulations
title: Simulate Sibling Pairs
sidebar_label: Sibling Pairs
---

:::note Requires PopGenSims.jl
To perfom simulations, you will need add and import the package `PopGenSims.jl` (available [here](https://github.com/pdimens/PopGenSims.jl)).
:::

It's good practice to use your data to simulate sibling pairs and perform
relatedness analyses on the simulations to understand the validity/strength 
of relatedness estimators on your data. To do this, you can use `simulate_sibship`
and specify the relationship you want to simulate and how many pairs to create for
that relationship.

```julia
simulate_sibship(data::PopData; n::Int, relationship::String, ploidy::Int)
```

This function will simulate mating crosses to generate `n` sample pairs (default: `500`) 
having the specified `relationship`, returning a `PopData` object. The simulations will first 
generate parents of the given `ploidy` (default: `2`) by drawing alleles from a global 
allele pool derived from the given `data` (i.e. weighted by their frequencies), then simulate mating between them.

### Relationship
Simulated parents will be crossed to generate offspring depending on the `relationship`:
- `"fullsib"` : 2 parents generate 2 full-sibling offspring, returns 2 offspring
- `"halfsib` : 3 parents generate 2 half-sibling offspring, returns 2 offspring
- `"unrelated"` : returns 2 randomly generated individuals from the global allele pools
- `"parent-offspring"` : 2 parents generate 1 offspring, returns 1 offspring and 1 parent

### Identifying pairs
The relationship between the newly generated samples can be identified by:
- Sample `name` will specify the simulation number, relationship, and whether parent or offspring
- Their `population` name will be that of their relationship (e.g. "fullsib")

:::tip plugging into relatedness
The `relatedness` function will recognize the population names output from simulating siblingship
and only estimate relatedness for the appropriate pairs. If you need this functionality, you are
strongly discouraged from manually editing the resulting `PopData` from `simulate_sibship`.
:::

### Ploidy
By default, the ploidy of the simulated parents and offspring are inferred from the supplied `PopData`.

:::note adjusting ploidy
If you have mixed-ploidy data or wish to generate parents and offspring of a ploidy different than the source
`PopData` you can specify the ploidy with which to simulate parents and offspring. For example, if your `PopData`
is diploid, but you wish to generate triploid or octoploid parents and offspring, you would specify `ploidy = 3`
 or `ploidy = 8` repectively.
:::

#### Odd ploidy
If trying to create offspring with an odd ploidy (3, 5, etc.), each parent has a 50% chance of 
contributing (½ × ploidy) + 1 alleles for all loci to the offspring. In other words, if ploidy = 3,
there's a 50% chance parent_1 will give 2 alleles for every locus for that simulated cross.

### Example
```
julia> cats = @nancycats ;

julia> fullsib_sims = simulate_sibship(cats, n = 50, relationship= "fullsib")
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent

julia> fullsib_sims.meta_df100×5 DataFrame
│ Row │ name          │ population │ ploidy │ longitude │ latitude │
│     │ String        │ String     │ Int64  │ Float32?  │ Float32? │
├─────┼───────────────┼────────────┼────────┼───────────┼──────────┤
│ 1   │ sim1_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
│ 2   │ sim1_fs_off2  │ fullsib    │ 2      │ missing   │ missing  │
│ 3   │ sim2_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
│ 4   │ sim2_fs_off2  │ fullsib    │ 2      │ missing   │ missing  │
│ 5   │ sim3_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
⋮
│ 95  │ sim48_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 96  │ sim48_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
│ 97  │ sim49_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 98  │ sim49_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
│ 99  │ sim50_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 100 │ sim50_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
```