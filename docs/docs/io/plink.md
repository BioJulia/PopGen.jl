---
id: plink
title: Plink
sidebar_label: Plink
---

## Import PLINK files as `PopData`
```julia
plink(infile::String; keepfields::Symbol|Vector{Symbol}, silent::Bool)
```
Read a PLINK `.ped` or binary `.bed` file into memory as a `PopData` object.
Requires an accompanying `.fam` file in the same directory, but an accompanying `.bim` file is optional.

- `infile::String` : path to `.ped` or `.bed` file

### Keyword Arguments
- `famfields::Symbol|Vector{Symbol}`: which additional fields to import from the `.fam` file
    - `:all` (default)
    - `:none`
    - any one or combination of `[:sire, :dam, :sex, :phenotype]`
- `bimfields::Symbol|Vector{Symbol}`: which additional fields to import from the optional `.bim` file
    - `:all` (default)
    - `:none`
    - any one or combination of `[:chromosome, :cm, :bp]`
- `silent::Bool`: whether to print file information during import (default: `false`)

**Example**

```julia
# assumes there is parakeet.ped + parakeet.fam in same directory
julia> parakeet = plink("datadir/parakeet.ped", famfields = :sex)

# assumes there is parrot.ped + parrot.fam in same directory
julia> parrot = plink("datadir/parrot.bed", famfields = [:sire, :dam], bimfields = :chromosome)
```

## Write PopData to PLINK format
```julia
plink(data::PopData; filename::String)
```
Write a biallelic `PopData` object to PLINK `.ped` format with an accompanying
`.fam` file. Genotypes are coded by the PLINK standard:
- Integers are the alleles
- `0` encodes missing
- After column 6, every two numbers indicate a diploid genotype such that:
    - `00` Homozygous for first allele
    - `01` Missing genotype
    - `10` Heterozygous
    - `11` Homozygous for second allele

**Example**

```julia
julia> sharks = dropmultiallelic(@gulfsharks) ;
julia> plink(sharks, filename = "biallelic_sharks.ped")
```