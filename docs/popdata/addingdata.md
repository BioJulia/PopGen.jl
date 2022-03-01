---
id: addingdata
title: Adding sample/locus data
sidebar_label: Adding sample/locus data
---

PopData's standard/default format includes information on:
- sample name
- sample population name
- sample ploidy
- sample genotypes
- locus name
- locus physical location (basepairs) [optional]
- locus genetic position (centiMorgans) [optional]

But, sometimes you might want to add more information to the data structure, so there are the mutating 
functions `sampleinfo!` and `locusinfo!` to do that.

## adding sampleinfo
```julia
sampleinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
sampleinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
```
Add an additional column sample information to `PopData` metadata. Mutates `PopData` in place. The new values 
must be in the same order as the samples in `sampleinfo(popdata)`.

### Arguments
- `metadata` : A Pair of `:ColumnName => [Values]`

### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)

**Example**
```
cats = @nancycats ;
sampleinfo!(cats, :whiskerlength => rand(metadata(cats).samples))
sampleinfo!(cats, "tailcolor" => rand(["orange", "brown"], metadata(cats).samples), categorical = true)
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["whiskerlength", "tailcolor"]
```
## adding locus information
```julia
locusinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
locusinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
```
Add an additional locus information to `PopData` metadata. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `locusinfo(PopData)`.

### Arguments
- `metadata` : A Pair of :ColumnName => [Values]

### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)

**Example**
```julia
cats = @nancycats
locusinfo!(cats, :quality => rand(metadata(cats).loci))
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["quality"]
```