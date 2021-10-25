---
id: ioutils
title: ioUtils.jl
sidebar_label: ioUtils.jl
---

### `determine_marker`
```julia
determine_marker(geno_parse::T, digits::Int) where T<:AbstractDataFrame
```
Return either `Int8` or `Int16` depending on largest allelic
value in all genotypes in the first 10 samples of an input 
DataFrame (or all the samples if less than 10 samples).
If the largest allele is 11 or greater, the marker will be 
considered a Microsatellite and coded in `PopData` as `Int16`, 
and the opposite is true for SNPs. There's no specific reason 
10 was chosen other than it being a reasonable buffer for edge
cases since SNP data <= 4, and haplotyped data could be a bit 
higher. Even if the microsatellite markers are coded 
incorrectly, there will be zero impact to performance,
and considering how few microsatellite markers are used in 
typical studies, the effect on in-memory size should be 
negligible (as compared to SNPs).

----

### `find_ploidy`
```julia
find_ploidy(genotypes::T where T<:SubArray)
```
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes of a sample and return the ploidy of the first non-missing locus.

----

### `phase`

```julia
phase(loc::Union{String, Int}, type::DataType, digit::Int)
```
Takes a String of numbers or Integer and returns a typed locus 
appropriate for PopGen.jl as used in the `genepop` and `csv` 
file parsers. Use `type` to specify output type (`Int8` or 
`Int16`), and `digit` to specify the number of digits/
characters used per allele in a locus.

**Example**

```julia
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
```
```julia
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```

----

### `unphase`

```julia    
unphase(geno::T; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) where T <: Genotype
unphase(geno::missing; digits::Int = 3, ploidy::Int = 2, miss::Int = 0) 
```
Takes a `Genotype` e.g. `(131, 94)` and returns a string of concatenated
alleles padded with *n* number of zeroes, where *n* is given by `digits =` .

#### miss
- `miss = 0`: `missing` values are returned as a string of `digits x ploidy` zeroes (default)
- `miss = -9` : `missing` values are returned as `-9`
The `ploidy` flag is only relevant for unphasing `missing` genotypes
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