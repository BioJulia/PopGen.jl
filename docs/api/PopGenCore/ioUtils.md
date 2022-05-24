---
id: ioutils
title: ioUtils.jl
sidebar_label: ioUtils.jl
---
## PopGenCore.jl/src/Utils/ioUtils.jl
| 📦  not exported | 🟪  exported by PopGenCore.jl | 🔵  exported by PopGen.jl |
|:---:|:---:|:---:|

### 📦 isbinary
```jula
isbinary(filepath::String)
```
Returns `true` if the `filepath` is a binary file. 

----

### 🟪 findploidy
```julia
findploidy(genotypes::T) where T<:AbstractVector
```
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes
of a sample and return the ploidy of the first non-missing locus.


----

### 🟪 phase
```julia
phase(loc::T, type::DataType, digit::Int) where T<:AbstractString
phase(loc::Missing, type::DataType, digit::Int) = missing
phase(loc::T, type::DataType, digits::T) where T<:Integer
```
Takes a String of numbers or Integers and returns a typed locus appropriate for PopGen.jl as used in the
`genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`),
and `digit` to specify the number of digits/characters used per allele in a locus.

**Examples**
```
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```

----

### 🟪 unphase
```julia
unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
unphase(geno::Missing; digits::Int = 3, ploidy::Int, miss::Int = 0)
```
Takes a `Genotype` (e.g. `(131, 94)`) and returns a string of concatenated
alleles padded with *n* number of zeroes, where *n* is given by `digits = `.
`missing` values are returned as either a string of 'digits × ploidy' zeroes (`miss = 0`)
or `"-9"` (`miss = -9`). The `ploidy` flag is only relevant for unphasing `missing` genotypes
and not used otherwise.

**Example**
```
unphase((1,2,3,4), digits = 3)
"001002003004"
unphase(missing, digits = 2, ploidy = 2, miss = -9)
"-9"
unphase(missing, digits = 2, ploidy = 2, miss = 0)
"0000"
```