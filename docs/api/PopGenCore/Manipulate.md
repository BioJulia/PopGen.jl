---
id: manipulate
title: Manipulate.jl
sidebar_label: Manipulate.jl
---
## PopGenCore.jl/src/Manipulate.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪🔵 sampleinfo!
```julia
sampleinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
sampleinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
```
Add an additional sample information to `PopData` metadata. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `sampleinfo(popdata)`.
#### Arguments
- `metadata` : A Pair of :ColumnName => [Values]
#### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)
**Example**
```
cats = @nancycats
sampleinfo!(cats, :whiskerlength => rand(cats.metadata.samples))
sampleinfo!(cats, "tailcolor" => rand(["orange", "brown"], metadata(cats).samples), categorical = true)
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["whiskerlength", "tailcolor"]
```

----
### 🟪🔵 locusinfo!
```julia
locusinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
locusinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
```
Add an additional locus information to `PopData` metadata. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `locusinfo(PopData)`.
#### Arguments
- `metadata` : A Pair of `:ColumnName => [Values]`
#### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)
**Example**
```
cats = @nancycats
locusinfo!(cats, :quality => rand(metadata(cats).loci))
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["quality"]
```
----
### 🟪🔵 locationdata!
```julia
locationdata!(data::PopData; longitude::Vector{Float64}, latitude::Vector{Float64})
```
Replaces existing `PopData` geographic coordinate data.
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
locationdata!(ncats, longitude = x, latitude = y)
```
----
```julia
locationdata!(data::PopData; longitude::Vector{String}, latitude::Vector{String})
```
Replaces existing `PopData` geographic coordinate data. Takes
**decimal minutes** or **degrees minutes seconds** format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
#### Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
##### NOTE
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
ncats = @nancycats;
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
locationdata!(ncats, longitude = x, latitude = y)
```

----
```julia
locationdata!(data::PopData, longitude::Vector{String}, latitude::Vector{String})
locationdata!(data::PopData; kwargs...)
```
----
### 🟪🔵 populations!
```
# Replace by matching
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Multiple methods to rename or reassign population names to a `PopData`.
#### Rename using a Dictionary
Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`

**Example**
```
potatopops = Dict("1" => "Idaho", "2" => "Russet")
populations!(potatoes, potatopops)
```
#### Rename using a Vector of Strings
If the number of new names is equal to the number of current unique population names,
the method will rename the existing populations in the order with which they appear 
via `unique()`. If the number of new population names is equal to the number of samples,
the method will instead assign new population names to every sample in the order with which they appear in `sampleinfo(popdata)`.

**Example**
```
# rename (2) existing populations
potatopops = ["Idaho", "Russet"]
populations!(potatoes, potatopops)
# assign new names to all [44] samples
potatopops = repeat(["Idaho", "Russet"], inner = 22) ;
populations!(potatoes, potatopops)
```
#### Reassign using samples and new population assignments
Completely reassign populations for each individual. Takes two vectors of strings
as input: one of the sample names, and the other with their new corresponding
population name. This can be useful to change population names for only some individuals.

**Example**
```
populations!(potatoes, ["potato_1", "potato_2"], ["north", "south"])
```
----
### 🟪🔵 exclude!
```julia
exclude!(data::PopData, kwargs...)
```
Edit a `PopData` object in-place by excluding all occurences of the specified information.
The keywords can be used in any combination. Synonymous with `omit!` and `remove!`. All
values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.
#### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
The keyword `loci` also works.
#### `population`
A `String` or `Vector{String}` of populations you want to remove from the `PopData`
The keyword `populations` also works.
#### `name`
A `String` or `Vector{String}` of samples you want to remove from the `PopData`
The keywords `names`, `sample`, and `samples` also work.

**Examples**
```
cats = @nancycats;
exclude!(cats, name = "N100", population = 1:5)
exclude!(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude!(cats, name = "N102", locus = :fca8, population = "3")
```
```julia
const omit! = exclude!
const remove! = exclude!
```
----
### 🟪🔵 exclude
```julia
exclude(data::PopData, kwargs...)
```
Returns a new `PopData` object excluding all occurrences of the specified keywords.
The keywords can be used in any combination. Synonymous with `omit` and `remove`. All
values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.
#### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
#### `population`
A `String` or `Vector{String}` of populations you want to remove from the `PopData`.
#### `name`
A `String` or `Vector{String}` of samples you want to remove from the `PopData`.

**Examples**
```
cats = @nancycats;
exclude(cats, name = "N100", population = 1:5)
exclude(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude(cats, name = "N102", locus = :fca8, population = "3")
```

```julia
const omit = exclude
const remove = exclude
```
----
### 🟪🔵 keep!
```julia
keep!(data::PopData, kwargs...)
```
Edit a `PopData` object in-place by keeping only the occurrences of the specified keywords.
If using multiple fields, they will be chained together as "`or`" statements.
All values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.
#### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to keep in the `PopData`.
#### `population`
A `String` or `Vector{String}` of populations you want to keep in the `PopData`.
#### `name`
A `String` or `Vector{String}` of samples you want to keep in the `PopData`.

**Examples**
```
cats = @nancycats;
keep!(cats, population = 1:5)
# keep 4 populations and 3 specific samples
keep!(cats, name = ["N100", "N102", "N211"])
# keep 2 loci, 2 populations, and 10 specific individuals
keep!(cats, locus = [:fca8, "fca37"], population = [7,8], name = samplenames(cats)[1:10])
```
----
### 🟪🔵 keep
```julia
keep(data::PopData, kwargs...)
```
Returns a new `PopData` object keeping only the occurrences of the specified keyword.
Unlike `exclude()`. only one keyword can be used at a time. All values are 
converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.
#### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to keep in the `PopData`.
#### `population`
A `String` or `Vector{String}` of populations you want to keep in the `PopData`.
#### `name`
A `String` or `Vector{String}` of samples you want to keep in the `PopData`.

**Examples**
```
cats = @nancycats;
keep(cats, population = 1:5)
# equivalent to cats[cats.genodata.population .∈ Ref(string.(1:5)), :]
keep(cats, name = ["N100", "N102", "N211"])
# equivalent to cats[cats.genodata.name .∈ Ref(["N100", "N102", "N211"]), :]
keep(cats, locus = [:fca8, "fca37"])
# equivalent to cats[cats.genodata.locus .∈ Ref(["fca8", "fca37"]), :]
```
----
### 🟪🔵 filter
```julia
filter(data::PopData, args...)
```
A drop-in replacement for DataFrames.filter where `PopData` is the first
argument and the filtering conditions are the second argument. Returns a
new `PopData`. **Note** the argument order is opposite of that from DataFrames.jl.
**Example**
```
x = @nancycats ;
y = filter(x, :name => i -> i ∈ samplenames(x)[1:10]) ;
show(x)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
show(y)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 10
  Populations: 1
```

----
### 🟪🔵 filter!
```julia
filter!(data::PopData, args...)
```
A drop-in replacement for the DataFrames.filter! where `PopData` is the first
argument and the filtering conditions are the second argument. Mutates the 
`PopData` in place and returns it. **Note** the argument order is opposite of that from DataFrames.jl.

**Example**
```
x = @nancycats ;
filter!(x, :name => i -> i ∈ samplenames(x)[1:10]) ;
show(x)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 10
  Populations: 1
```
