## Import a genepop file as a `PopObj`

```julia
genepop(infile; digits = 3, popsep = "POP", marker = "snp", diploid = true)

# Example
julia> b = genepop("/data/wasp_hive.gen", digits = 3, popsep = "POP")
```

??? warning "Windows users"
    make sure to change your backslashes "\" to forward slashes "/" 

### arguments

- `#!julila infile::String` : path to genepop file, in quotes

### keyword Arguments

- `#!julila digits::Int64` : the number of digits used to denote an allele (default: `3`)
- `#!julila popsep::String` : word that separates populations in `infile` (default: `"POP"`)
- `#!julila marker::String` : "snp" (default) or "msat" for microsatellites
- `#!julila diploid::Bool` :  uses memory-optimized parsing for diploid samples (default: `true`)

!!! info ""
    By default, the file reader will assign numbers as population ID's (as Strings) in order of appearance in the genepop file. Use the `populations!` function to rename these with your own population ID's.

## Format

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is **tab** delimited

### formatting examples

```tab="loci stacked vertically"
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

```tab="loci stacked horizontally"
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



## Acknowledgements

The original implementations of this parser were written using only Base Julia, and while the speed was fantastic, the memory footprint involved seemed unusually high (~650mb RAM to parse `gulfsharks`, which is 3.2mb in size). However, thanks to the efforts of the [CSV.jl](https://github.com/JuliaData/CSV.jl) and [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) teams, we leverage those packages to do much of the heavy lifting, while preserving the speed and reducing the memory footprint by as much as 65% in diploid cases (~166mb RAM to process `gulfsharks`). 

