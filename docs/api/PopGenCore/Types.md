---
id: types
title: Types.jl
sidebar_label: Types.jl
---

### `PopObj`
```Julia
PopObj
```
Generic AbstractType for use in PopGen.jl

----

### `PopData`
```
PopData(meta::IndexedTable, loci::IndexedTable)
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into PopData,
and instead should use the provided genepop, csv, or vcf file importers.

**`meta` ::DataFrame** individual/sample data with the columns:

- `name` ::String the individual/sample names
- `population` ::String population names
- `ploidy` ::Int8 ploidy in order of `ind`
- `longitude` ::Float64 longitude values
- `latitude` ::Float64 latitude values

**`loci` ::DataFrame** Long-format table of sample genotype records

- `name` ::String of the individual/sample names
- `population`::String population name
- `locus` ::String of locus name
- `genotype` Tuple of Int8 or Int16 depending on SNP or microsatellite

----

### `GenoType`
```julia
Genotype::DataType
```
For convenience purposes, an alias for `NTuple{N, <:Signed} where N`, which is the type describing individual genotypes in PopData.

----

### `GenoTypeArray`
```julia
GenotypeArray::DataType
```
For convenience purposes, an alias for an `AbstractVector` of elements `Missing` and `Genotype`, which itself is of type `NTuple{N, <:Signed} where N`. The definition as an `AbstractVector` adds flexibility for `SubArray` cases.

----

### `show`
    Base.show(io::IO, data::PopData)
Overloads `Base.show` for concise summary printing of a PopData object.