---
id: manipulate
title: Manipulate.jl
sidebar_label: Manipulate.jl
---
## PopGenCore.jl/src/Manipulate.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### `add_meta!`

```julia
add_meta!(popdata::PopData, metadata::T; name::String, loci::Bool = true, categorical::Bool = true) where T <: AbstractVector
```
Add an additional metadata information to a `PopData` object. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `PopData.sampleinfo`.

#### Arguments
- `popdata` : The `PopData` object to add information to
- `metadata` : A `Vector` with the metadata you wish to add to the `PopData`, in the same order as the names appear in `PopData.sampleinfo`

#### Keyword Arguments
- `name` : String of the name of this new column
- `loci` : Boolean of whether to also add this information to `PopData.genodata` (default: `true`)
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `true`)

----

### `add_meta!`
```julia
add_meta!(popdata::PopData, samples::Vector{String}, metadata::T; name::String, loci::Bool = true, categorical::Bool = true) where T <: AbstractVector
```
Add an additional metadata information to a `PopData` object. Mutates `PopData` in place. Takes a vector of
sample names if the metadata is not in the same order as samples appear in `PopData.sampleinfo`.

#### Arguments
- `popdata` : The `PopData` object to add information to
- `sample` : A `Vector{String}` of sample names corresponding to the order of the `metadata` 
- `metadata` : A `Vector` with the metadata you wish to add to the `PopData`, in the same order as the names appear in `PopData.sampleinfo`

#### Keyword Arguments
- `name` : String of the name of this new column
- `loci` : Boolean of whether to also add this information to `PopData.genodata` (default: `true`)
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `true`)

----
### `locations`
```julia
locationdata(data::PopData)
```
View the longitude and lattitude data in a `PopData` object. Returns a table derived from the PopData. Changes made to this table will not alter the source `PopData` object.

Use `locations!` to add spatial data to a `PopData` object.

----

### `locations!`
```julia
locationdata!(data::PopData; longitude::Vector{Float64}, lattitude::Vector{Float64})
```
Replaces existing `PopData` location data (longitude `long`, lattitude `lat`).
Takes **decimal degrees** as a `Vector` of any `AbstractFloat`.

#### Formatting requirements
- Decimal Degrees format: `-11.431`
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)

**Example**
```
ncats = @nancycats ;
x = rand(237) ; y = rand(237)
locationdata!(ncats, longitude = x, lattitude = y)
```

----

### `locations!`
```julia
locationdata!(data::PopData; longitude::Vector{String}, lattitude::Vector{String})
```
Replaces existing `PopData` location data (longitude `long`, lattitude `lat`). Takes
**decimal minutes** format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.

#### Formatting requirements
- Decimal Minutes: `"-11 43.11"` (must use space and be a `String`)
- **Must** use negative sign `-` or single-letter cardinal directions like "11 43.11W"
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)

**NOTE**

If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, lattitude degrees, lattitude minutes),
then the easiest course of action would be to merge them into two vectors of strings
(one for longitude, one for lattitude):
```
long_string = string.(lat_deg, " ", lat_min)
lat_string = string.(long_deg, " ", long_min)
```
and use these as inputs into `locations!`

**Example**
```
ncats = @nancycats;
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
locationdata!(ncats, longitude = x, lattitude = y)
```


----

### `loci`
```julia
loci(data::PopData)
```
Returns an array of strings of the loci names in a `PopData` object.

----

### `get_genotypes`
```julia
get_genotypes(data::PopObj; samples::Union{String, Vector{String}}, loci::Union{String, Vector{String}})
```
Return the genotype(s) of one or more `samples` for one or more specific `loci` (both as keywords) in a `PopData` object.

**Example**
```julia
cats = @nancycats;
get_genotype(cats, samples = "N115" , loci = "fca8")
get_genotypes(cats, samples = ["N1", "N2"] , loci = "fca8")
get_genotype(cats, samples = "N115" , loci = ["fca8", "fca37"])
get_genotype(cats, samples = ["N1", "N2"] , loci = ["fca8", "fca37"])
```

----

### `get_sample_genotypes`
```julia
get_sample_genotypes(data::PopData, sample::String)
```
Return all the genotypes of a specific sample in a `PopData` object. This is an extension for the internal function `get_genotypes`.

**Example**
```julia
cats = @nancycats
get_sample_genotypes(cats, "N115")
```

----

### `locus`
```julia
locus(data::PopData, locus::Union{String, Symbol})
```
Convenience wrapper to return a vector of all the genotypes of a single locus

**Example**
```julia
locus(@gulfsharks, "contig_475")
```

----

### `populations`

```julia
    populations(data::PopData; counts::Bool = false)
```
View unique population ID's in a `PopData` objec

- `counts = true` returns the number of samples per population  (default = `false`)

----

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
`Vector` of new unique population names in the order that they appear in the PopData.sampleinfo

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

----

### `exclude!`
```julia
exclude!(data::PopData, kwargs...)
remove!(data::PopData, kwargs...)
omit!(data::PopData, kwargs...)
```
Edit a `PopData` object in-place by excluding all occurences of the specified information.
The keywords can be used in any combination. Synonymous with `omit!` and `remove!`.

**Keyword Arguments**

-`locus`: A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
    - The keyword `loci` also works.
- `population`: A `String` or `Vector{String}` of populations you want to remove from the `PopData`.
    - The keyword `populations` also works.
- `name`: A `String` or `Vector{String}` of samples you want to remove from the `PopData`.
    - The keywords `names`, `sample`, and `samples` also work.

**Examples**

```julia
cats = @nancycats;
exclude!(cats, name = "N100", population = ["1", "15"])
exclude!(cats, samples = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude!(cats, names = "N102", loci = "fca8", population = "3")
```

----

### `exclude`
```julia
exclude(data::PopData, kwargs...)
remove(data::PopData, kwargs...)
omit(data::PopData, kwargs...)
```
Returns a new `PopData` object excluding all occurrences of the specified keywords.
The keywords can be used in any combination. Synonymous with `omit` and `remove`.

**Keyword Arguments**

- `locus`: A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
    - The keyword `loci` also works.
- `population`: A `String` or `Vector{String}` of populations you want to remove from the `PopData`.
    - The keyword `populations` also works.
- `name`: A `String` or `Vector{String}` of samples you want to remove from the `PopData`.
    - The keywords `names`, `sample`, and `samples` also work.

**Examples**

```julia
cats = @nancycats;
exclude(cats, name = "N100", population = ["1", "15"])
exclude(cats, samples = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude(cats, names = "N102", loci = "fca8", population = "3")
```
