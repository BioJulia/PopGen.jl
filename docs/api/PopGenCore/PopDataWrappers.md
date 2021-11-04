---
id: popdatawrappers
title: PopDataWrappers.jl
sidebar_label: PopDataWrappers.jl
---
📦  => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 🟪🔵 info
```julia
info(::PopData)
```
Display the metadata (PopDataInfo) of a PopData object.

----
### 🟪🔵 metadata
```julia
metadata(::PopData)
```
Display the metadata (PopDataInfo) of a PopData object.

----
### 🟪🔵 sampleinfo
```julia
sampleinfo(::PopData)
```
Show the sample information found within the metadata of a `PopData` object. Returns a view of a dataframe

----
### 🟪🔵 genodata
```julia
genodata(::PopData)
```
Method to show the genotype information of a `PopData` object. Returns a view of a dataframe. 

----
### 🟪🔵 locusinfo
```julia
locusinfo(::PopData)
```
Show the locus information found within the metadata of a `PopData` object. Returns a view of a dataframe

----
### 🟪🔵 locationdata
```julia
locationdata(data::PopData)
```
View the longitude and latitude data in a `PopData` object. Returns a table
derived from the PopData. Changes made to this table will not alter the source
`PopData` object.
Use `locations!` to add spatial data to a `PopData` object.

----
### 🟪🔵 loci
```julia
loci(data::PopData)
```
Returns an array of strings of the loci names in a `PopData` object.

----
### 🟪🔵 samplenames
```julia
samplenames(data::PopData)
```
View individual/sample names in a `PopData`

----
### 🟪🔵 genotypes
```julia
genotypes(data::PopData, samplelocus::String)
```
Return a vector of all the genotypes of a sample (or locus) in a `PopData` object.
```
cats = @nancycats
genotypes(cats, "N115")
genotypes(cats, "fca8")
```

----
```julia
genotypes(data::PopData, samplelocus::Pair{String, String}) ::DataFrame
genotypes(data::PopData, samplelocus::Pair{Vector{String}, String}) ::DataFrame
genotypes(data::PopData, samplelocus::Pair{String, Vector{String}}) ::DataFrame
```    
Return a genotype or dataframe of genotypes for one or more samples/loci 
in a `PopData` object. Uses the `Pair` notation of `samples => loci`.

**Examples**
```
cats = @nancycats;
genotypes(cats, "N115" => "fca8")
genotypes(cats, ["N115", "N7"] => "fca8")
genotypes(cats, "N115" => ["fca8", "fca37"])
genotypes(cats, ["N1", "N2"] => ["fca8", "fca37"])
```

----
### 🟪🔵 populations
```julia
populations(data::PopData; counts::Bool = false)
```
View unique population ID's and/or their counts in `PopData`.
- `counts` returns a dataframe of samples per `population` instead (default = `false`)

