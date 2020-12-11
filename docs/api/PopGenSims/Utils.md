---
id: popgensims_utils
title: Utils.jl
sidebar_label: Utils.jl
---

### `append!`
```julia
append!(data::PopData, data2::PopData)
```
Add the rows of `data2` to the end of `data`. This will add the samples present
in the second `PopData` object to the first `PopData` object (mutating it). 
**Note** that this is a simple appending, and you risk corrupting your `PopData` if
the two `PopData` objects do not have identical loci.

**Example**
```
julia> cats = @nancycats
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 237
  Loci: 9
  Populations: 17
  Longitude: absent
  Latitude: absent

julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent

julia> append!(cats, purrfect_pairs)
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 337
  Loci: 9
  Populations: 18
  Longitude: absent
  Latitude: absent
```
----
### `append`

```julia
append(data::PopData, data2::PopData)
```
Add the rows of `data2` to the end of `data`. This will combine the samples present
in both `PopData` objects and return a new `PopData` object. **Note** that this is 
a simple appending, and you risk corrupting your `PopData` if the two `PopData` 
objects do not have identical loci.

**Example**
```
julia> cats = @nancycats
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 237
  Loci: 9
  Populations: 17
  Longitude: absent
  Latitude: absent

julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent

julia> merged_cats = append(cats, purrfect_pairs)
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 337
  Loci: 9
  Populations: 18
  Longitude: absent
  Latitude: absent
```

-----

### `simulate_sample`
```julia
simulate_sample(alleles::Dict{String,NTuple}, loc::Vector{String}; ploidy::Int)
```
Using a global allele pool given by a Dict{loci,alleles} and a list of loci (`loc`), simulate
an individual with a given `ploidy`. Returns a Vector of genotypes.

**Example**
```
julia> cats = @nancycats ;
julia> loc, alleles = allele_pool(cats) ;
julia> simulate_parent(alleles, loc, ploidy = 2)
9-element Array{Array{Int16,1},1}:
 [139, 129]
 [146, 146]
 [145, 141]
 [126, 126]
 [150, 148]
 [148, 140]
 [185, 199]
 [91, 113]
 [208, 208]
```
