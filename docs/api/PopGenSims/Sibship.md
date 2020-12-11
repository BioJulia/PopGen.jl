---
id: popgensims_sibship
title: Sibship.jl
sidebar_label: Sibship.jl
---

```julia
_cross(parent1::Vector{Vector{T}}, parent2::Vector{Vector{T}}) where T <: Signed
```
Simulate a mating cross between two parents, generating one offspring with the same
ploidy as `parent1`. This variant of `cross` is used internally for `simulate_sibship`.

----

### `parentoffspring`
```julia
parentoffspring(data::PopData; n::Int, ploidy::Int)
```

----

### `fullsib`
```julia
fullsib(data::PopData; n::Int, ploidy::Int)
```
----

### `halfsib`
```julia
halfsib(data::PopData; n::Int, ploidy::Int)
```
----

### `unrelated`
```julia
unrelated(data::PopData; n::Int, ploidy::Int)
```

----

### `simulate_sibship`
```julia
simulate_sibship(data::PopData; n::Int, relationship::String, ploidy::Int)
```
Simulate mating crosses to generate `n` sample pairs (default: `500`) having the specified `relationship`, 
returning a `PopData` object. The simulations will first generate parents of the given `ploidy` (default: `2`) 
by drawing alleles from a global allele pool derived from the given `data` (i.e. weighted by their frequencies).

#### Relationship
Simulated parents will be crossed to generate offspring depending on the `relationship`:
- `"fullsib"` : 2 parents generate 2 full-sibling offspring, return 2 offspring
- `"halfsib` : 3 parents generate 2 half-sibling offspring, returns 2 offspring
- `"unrelated"` : returns 2 randomly generated individuals from the global allele pools
- `"parent-offspring"` : 2 parents generate 1 offspring, returns 1 offspring and 1 parent

#### Identifying pairs
The relationship between the newly generated samples can be identified by:
- Sample `name`s will specify their simulation number, relationship, and whether parent or offspring
- Their `population` name will be that of their relationship (e.g. "fullsib")

#### Ploidy
If your data is not diploid, then change this value to the appropriate ploidy. While the simulations default to 
diploid, if you wish to generate parents and offspring of a ploidy different than the source `PopData` you can 
change this value. For example, if your `PopData` is diploid, but you wish to generate triploid or octoploid 
parents and offspring, you can. 
#### Odd ploidy
If trying to create offspring with an odd ploidy (3,5, etc.), each parent has a 50% chance of 
contributing (½ × ploidy) + 1 alleles for all loci to the offspring. In other words, if ploidy = 3,
there's a 50% chance parent_1 will give 2 alleles for every locus for that simulated cross.

**Example**
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