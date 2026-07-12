---
id: plink
title: Plink.jl
sidebar_label: Plink.jl
---
## PopGenCore.jl/src/io/Plink.jl
| 📦  not exported | 🟪  exported by PopGenCore.jl | 🔵  exported by PopGen.jl |
|:---:|:---:|:---:|

### 📦 _plinkindex
```julia
_plinkindex(s::Matrix{UInt8}, i::Integer, j::Integer)
_plinkindex(s::Matrix{UInt8})
```
Tthis function copies the getindex tool OpenMendel/SnpArrays.jl uses 
to pull out the byte values from the compressed hex genotypes
Repo: https://github.com/OpenMendel/SnpArrays.jl

----
### 📦 _SNP
```julia
_SNP(genotype::UInt8)
_SNP(genomatrix::AbstractArray{UInt8})
```
- 00	Homozygous for first allele (0x00)
- 01	Missing genotype (0x01)
- 10	Heterozygous  (0x02)
- 11	Homozygous for second allele in .bim file (0x03)
----
### 📦 _plinkped
```julia
_plinkped(infile::String, keepfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
```
----
### 📦 _plinkbed
```julia
_plinkbed(infile::String, famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
```
----
### 📦 _genoconversion
```julia
_genoconversion(genotype::T) where T<:Genotype = join(genotype, " ")
_genoconversion(genotype::Missing) = "0 0"
```

----
### 🟪🔵 plink
```julia
plink(infile::String; famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
```
Read a PLINK `.ped` or binary `.bed` file into memory as a `PopData` object.
Requires an accompanying `.fam` file in the same directory, but an accompanying `.bim` file is optional.
- `infile::String` : path to `.ped` or `.bed` file
### Keyword Arguments
- `famfields::Symbol|Vector{Symbol}` : which additional fields to import from the `.fam` file
    - `:all` [default]
    - `:none`
    - any one or combination of `[:sire, :dam, :sex, :phenotype]`
- `bimfields::Symbol|Vector{Symbol}` : which additional fields to import from the optional `.bim` file
    - `:all` [default]
    - `:none`
    - any one or combination of `[:chromosome, :cm, :bp]`
- `silent::Bool`   : whether to print file information during import (default: `false`)
## Example
```julia
parakeet = plink("datadir/parakeet.ped", famfields = :sex)
parrot = plink("datadir/parrot.bed", famfields = [:sire, :dam], bimfields = :chromosome)
```

----
```julia
plink(data::PopData; filename::String)
```
Write a biallelic `PopData` object to PLINK `.ped` format with an accompanying
`.fam` file. Genotypes are coded by the PLINK standard:
- Integers are the alleles
- `0` encodes missing
- After column 6, every two numbers indicate a diploid genotype.
## Example
```julia
sharks = dropmultiallelic(@gulfsharks) ;
plink(sharks, filename = "biallelic_sharks.ped")
```