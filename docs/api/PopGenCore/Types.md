---
id: types
title: PopData.jl
sidebar_label: PopData.jl
---
## PopGenCore.jl/src/PopData.jl
| 📦  not exported | 🟪  exported by PopGenCore.jl | 🔵  exported by PopGen.jl |
|:---:|:---:|:---:|

### 🟪🔵 PopObj
```Julia
PopObj::AbstractType
```
Generic AbstractType for use in PopGen.jl

### 🟪🔵 PopDataInfo
```
mutable struct PopDataInfo
    samples::Int64
    sampeinfo::DataFrame
    loci::Int64
    locusinfo::DataFrame
    populations::Int64
    ploidy::Union{Int8, Vector{Int8}}
    biallelic::Bool
end
```
The data struct used internally as `PopData.metadata` fields to store basic information
about the `PopData` for easy access.

```julia
PopDataInfo(genodf::DataFrame)
```
constructor format using just the genodata dataframe

### 🟪🔵 PopData
```
PopData
    metadata::PopDataInfo
    genodata::DataFrame
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into a PopData,
and instead should use the provided file importers and utilities.
    - `metadata` PopDataInfo of  data information
        - `samples` - the number of samples in the data
        - `sampleinfo` - DataFrame of sample names,populations, ploidy, etc.
        - `loci` - the number of loci in the data
        - `locusinfo` - DataFrame of locus names, chromosome, physical position, etc.
        - `populations` - the number of populations in the data
        - `ploidy` - the ploidy (or ploidies) present in the data
        - `biallelic` - if all the markers are biallelic
    - `genodata` DataFrame of sample genotype records
        - `name` - the individual/sample names [`PooledArray`]
        - `population` - population names [`PooledArray`]
        - `locus` - locus names [`PooledArray`]
        - `genotype` - genotype values [`NTuple{N,Signed}`]

```PopData(data::DataFrame)```
Contructor using just a `genodata` dataframe.

### 🟪🔵 PopDataInfo!
```julia
PopDataInfo!(data::PopData)
```
This method is used to update PopDataInfo from PopData, all in one swoop

```julia
PopDataInfo!(popdatainfo::PopDataInfo, genodata::DataFrame)
```
Method to update preexisting PopDataInfo with new genodata, useful for getindex and creating new PopData from that

### 🟪🔵 Genotype
```julia
Genotype::DataType
```
For convenience purposes, an alias for `NTuple{N, <:Signed} where N`, which is
the type describing individual genotypes in PopData. Specifically, there exist
`SNP` as an alias for `NTuple{N, Int8}` and `MSat` for `NTuple{N, Int16}`

### 🟪🔵 SNP
```julia
SNP::DataType
```
An alias for `NTuple{N, Int8}`

### 📦 _SNP
```julia
SNP(geno)
```
Contstructor for `SNP`

### 🟪🔵 MSat
```julia
MSat::DataType
```
An alias for `NTuple{N, Int16}`

### 📦 _MSat
```julia
_MSat(geno)
```
Constructor for `MSat`

### 🟪🔵 GenoArray
```julia
GenoArray::DataType
```
An alias for an `AbstractVector` of elements `Missing`
and `Genotype`, which itself is of type `NTuple{N, <:Integer} where N`.
The definition as an `AbstractVector` adds flexibility for `SubArray`
cases.

### 📦 _ploidy2text
```
_ploidy2text(ploidy::Int8)
_ploidy2text(ploidy::Vector{Int8})
```

### 🟪🔵 Base.show
```julia
Base.show(io::IO, data::PopData)
Base.show(io::IO, data::PopDataInfo)
```

### 🟪🔵 Base.getindex
```julia
Base.getindex(data::PopData, idx::Symbol)
Base.getindex(data::PopData, args)
Base.getindex(data::PopData, expression, cols)
```

### 🟪🔵 Base.getproperty
```julia
getproperty(data::PopData, field::Symbol)
```
A convenience method to access certain elements in a `PopData` with fewer keystrokes. 
Essentially a standard `getproperty` call, except `sampleinfo` accesses `metadata.sampleinfo`,
`locusinfo` accesses `metadata.locusinfo`, and `info` is an alias for `metadata`.

**Example**
```julia
cats = @nancycats ;
cats.metadata == cats.info
cats.metadata.sampleinfo == cats.sampleinfo
cats.metadata.locusinfo == cats.locusinfo
```
