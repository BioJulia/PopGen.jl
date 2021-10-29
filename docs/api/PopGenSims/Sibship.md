---
id: popgensims_sibship
title: Sibship.jl
sidebar_label: Sibship.jl
---
PopGenSims.jl/src/Sibship.jl
❗ => not exported | 
⚫ => exported by PopGenSims.jl


### ❗_cross
```julia
_cross(parent1::Vector{Vector{T}}, parent2::Vector{Vector{T}}) where T <: Signed
```
Simulate a mating cross between two parents, generating one offspring with the same
ploidy as `parent1`. This variant of `cross` is used internally for `simulate_sibship`.

----

### ❗_parentoffspring
```julia
parentoffspring(data::PopData; n::Int, ploidy::Int)
```

----

### ❗fullsib
```julia
fullsib(data::PopData; n::Int, ploidy::Int)
```
----

### ❗halfsib
```julia
halfsib(data::PopData; n::Int, ploidy::Int)
```
----

### ❗unrelated
```julia
unrelated(data::PopData; n::Int, ploidy::Int)
```

----

### ⚫ simulate_sibship
```julia
simulate_sibship(data::PopData; fullsib::Int, halfsib::Int, unrelated::Int, parentoffspring::Int, ploidy::Signed)
```
Simulate mating crosses to generate sample pairs with any combination of the specified relationships, 
returning a `PopData` object. The simulations will first generate parents of a given
`ploidy` (inferred or specified) by drawing alleles from a global allele pool derived
from the given `data` (i.e. weighted by their frequencies).

#### Relationship
Simulated parents will be crossed to generate offspring depending on the relationship:
- `fullsib` : 2 parents generate 2 full-sibling offspring, return 2 offspring
- `halfsib` : 3 parents generate 2 half-sibling offspring, returns 2 offspring
- `unrelated` : returns 2 randomly generated individuals from the global allele pools
- `parentoffspring` : 2 parents generate 1 offspring, returns 1 offspring and 1 parent

**Identifying pairs**

The relationship between the newly generated samples can be identified by:
- Sample `name`s will specify their simulation number, relationship, and whether parent or offspring
    - Naming convention: [simulation #]_[relationship]_[offspring #]
    - example: sim005_fullsib_1 = [simulation 005]_[full sibling]_[offspring 1]
- Their `population` name will be that of their relationship (e.g. "fullsib")

**Ploidy**

If the samples in your `PopData` are of a single ploidy, then `ploidy = 0` (the default) will infer the ploidy
and generate parents and offspring according to the ploidy of your data. If you have mixed-ploidy data or wish 
to generate parents and offspring of a ploidy different than the source `PopData` you can specify the ploidy
with which to simulate parents and offspring. For example, if your `PopData` is diploid, but you wish to generate
triploid or octoploid parents and offspring, you would specify `ploidy = 3` or `ploidy = 8` repectively. 

**Odd ploidy**

If trying to create offspring with an odd ploidy (3,5, etc.), each parent has a 50% chance of 
contributing (½ × ploidy) + 1 alleles for all loci to the offspring. In other words, if ploidy = 3,
there's a 50% chance parent_1 will give 2 alleles for every locus for that simulated cross.

**Example**
```
julia> cats = @nanycats ;

julia> cat_sims = simulate_sibship(cats, fullsib = 10, halfsib = 50)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 120
  Populations: 2

julia> cat_sims.sampleinfo
120×3 DataFrame
 Row │ name             population  ploidy 
     │ String           String      Int64  
─────┼─────────────────────────────────────
   1 │ sim01_fullsib_1  fullsib          2
   2 │ sim01_fullsib_2  fullsib          2
   3 │ sim02_fullsib_1  fullsib          2
   4 │ sim02_fullsib_2  fullsib          2
   5 │ sim03_fullsib_1  fullsib          2
   6 │ sim03_fullsib_2  fullsib          2
  ⋮  │        ⋮             ⋮         ⋮
 115 │ sim48_halfsib_1  halfsib          2
 116 │ sim48_halfsib_2  halfsib          2
 117 │ sim49_halfsib_1  halfsib          2
 118 │ sim49_halfsib_2  halfsib          2
 119 │ sim50_halfsib_1  halfsib          2
 120 │ sim50_halfsib_2  halfsib          2
                           108 rows omitted
```