---
id: genepop
title: Genepop format
sidebar_label: Genepop format
---

## Import a genepop file as `PopData`

```julia
genepop(infile; kwargs...)

# Example
julia> b = genepop("/data/wasp_hive.gen", digits = 3, popsep = "POP")
```

:::caution Windows users
Make sure to change your backslashes `\` to forward slashes `/` 
:::

### arguments

- `infile::String` : path to genepop file, in quotes

### keyword Arguments

- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent::Bool`   : whether to print file information during import (default: `true`)

:::info population names
By default, the file reader will assign numbers as population ID's (as Strings) in order of appearance in the genepop file. Use the `populations!` function to rename these with your own population ID's.
:::

## Format

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular and consistent keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is **tab** delimited or **space** delimited, but not both

:::: tabs card stretch
::: tab genepop w/loci stacked vertically
```
Wasp populations in New York
Locus1
Locus2
Locus3
POP
Oneida_01,	250230	564568	110100
Oneida_02,	252238	568558	100120
Oneida_03,	254230	564558	090100
POP
Newcomb_01,	254230	564558	080100
Newcomb_02,	000230	564558	090080
Newcomb_03,	254230	000000	090100
Newcomb_04,	254230	564000	090120
```
:::
::: tab genepop w/loci stacked horizontally
```
Wasp populations in New York
Locus1,Locus2,Locus3
POP
Oneida_01,	250230	564568	110100
Oneida_02,	252238	568558	100120
Oneida_03,	254230	564558	090100
POP
Newcomb_01,	254230	564558	080100
Newcomb_02,	000230	564558	090080
Newcomb_03,	254230	000000	090100
Newcomb_04,	254230	564000	090120
```
:::
::::

## Acknowledgements

The original implementations of this parser were written using only Base Julia, and while the speed was fantastic, the memory footprint involved seemed unusually high (~650mb RAM to parse `gulfsharks`, which is only 3.2mb in size). However, thanks to the efforts of [CSV.jl](https://github.com/JuliaData/CSV.jl) , we leverage that package to preserve the speed and reducie the memory footprint quite a bit!

