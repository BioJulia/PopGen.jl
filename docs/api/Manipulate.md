---
id: manipulate
title: Manipulate.jl
sidebar_label: Manipulate.jl
---

### `locations`
```julia
    locations(data::PopData)
```
View the longitude and latitude data in a `PopData` object. Returns a table derived from the PopData. Changes made to this table will not alter the source `PopData` object.

Use `locations!` to add spatial data to a `PopData` object.

### `locations!`
```julia
locations!(data::PopData; long::Vector{Float64}, lat::Vector{Float64})
```
Replaces existing `PopData` location data (longitude `long`, latitude `lat`).
Takes **decimal degrees** as a `Vector` of any `AbstractFloat`.
#### Formatting requirements
- Decimal Degrees format: `-11.431`
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)

**Example**
```
ncats = nancycats() ;
x = rand(237) ; y = rand(237)
locations!(ncats, long = x, lat = y)
```

### `locations!`
```julia
locations!(data::PopData; long::Vector{String}, lat::Vector{String})
```
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes
**decimal minutes** format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
#### Formatting requirements
- Decimal Minutes: `"-11 43.11"` (must use space and be a `String`)
- **Must** use negative sign `-` or single-letter cardinal directions like "11 43.11W"
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)

**NOTE**

If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, latitude degrees, latitude minutes),
then the easiest course of action would be to merge them into two vectors of strings
(one for longitude, one for latitude):
```
long_string = string.(lat_deg, " ", lat_min)
lat_string = string.(long_deg, " ", long_min)
```
and use these as inputs into `locations!`

**Example**
```
ncats = nancycats();
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
locations!(ncats, long = x, lat = y)
```


### `loci`
```julia
loci(data::PopData)
```
Returns an array of strings of the loci names in a `PopData` object.

### `get_genotypes`
```julia
get_genotypes(data::PopObj; samples::Union{String, Vector{String}}, loci::Union{String, Vector{String}})
```
Return the genotype(s) of one or more `samples` for one or more specific `loci` (both as keywords) in a `PopData` object.

**Example**
```julia
cats = nancycats();
get_genotype(cats, samples = "N115" , loci = "fca8")
get_genotypes(cats, samples = ["N1", "N2"] , loci = "fca8")
get_genotype(cats, samples = "N115" , loci = ["fca8", "fca37"])
get_genotype(cats, samples = ["N1", "N2"] , loci = ["fca8", "fca37"])
```

### `get_sample_genotypes`
```julia
get_sample_genotypes(data::PopData, sample::String)
```
Return all the genotypes of a specific sample in a `PopData` object. This is an extension for the internal function `get_genotypes`.

**Example**
```julia
cats = nancycats()
get_sample_genotypes(cats, "N115")
```

### `locus`
```julia
locus(data::PopData, locus::Union{String, Symbol})
```
Convenience wrapper to return a vector of all the genotypes of a single locus 

**Example**
```julia
locus(gulfsharks(), "contig_475")
```

### `missing`
```julia
missing(data::PopData; by::String = "sample")
```
Get missing genotype information in a `PopData`. Specify a mode of operation to return a DataFrame corresponding with that missing information.

**Modes**
- `"sample"` - returns a count and list of missing loci per individual (default)
- `"pop"` - returns a count of missing genotypes per population
- `"locus"` - returns a count of missing genotypes per locus
- `"full"` - returns a count of missing genotypes per locus per population

**Example**
```julia
missing(gulfsharks(), by = "pop")
```

```julia
    populations(data::PopData; listall::Bool = false)
```
View unique population ID's and their counts in a `PopData`.

- `listall = true` displays all samples and their `population` instead (default = `false`)

### `populations!`
```julia
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Multiple methods to rename or reassign population names in`PopData`.
#### Rename using a Dictionary
```julia
populations!(data::PopData, rename::Dict)
```
Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`

**Example**
```
potatopops = Dict("1" => "Idaho", "2" => "Russet")
populations!(potatoes, potatopops)
```
#### Rename using a Vector of Strings
```julia
populations!(data::PopData, rename::Vector{String})
```
`Vector` of new unique population names in the order that they appear in the PopData.meta

**Example**
```
potatopops = ["Idaho", "Russet"]
populations!(potatoes, potatopops)
```
#### Reassign using samples and new population assignments
```julia
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Completely reassign populations for each individual. Takes two vectors of strings
as input: one of the sample names, and the other with their new corresponding
population name.

**Example**
```
populations!(potatoes, ["potato_1", "potato_2"], ["north_russet", "south_russet"])
```

### `exclude_loci`
```julia
exclude_loci(data::PopData, locus::String)
exclude_loci(data::PopData, loci::Vector{String})
```
Exclude selected loci from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_loci` and `remove_loci`.

**Example**
```julia
new_cats = exclude_loci(nancycats(), "fca8")
very_new_cats = exclude_loci(nancycats(), ["fca8", "fca23"])
```

### `exclude_samples`
```julia
exclude_samples(data::PopData, samp_id::String)
exclude_samples(data::PopData, samp_id::Vector{String})
```
Exclude selected samples from a `PopData` object. Returns a new `PopData` object, leaving the original intact. Synonymous with `omit_samples` and `remove_samples`.

**Example**
```julia
exclude_samples(nancycats, "N100")
exclude_samples(nancycats, ["N100", "N102", "N211"])
```

```julia
samples(data::PopData)
```
View individual/sample names in a `PopData`
