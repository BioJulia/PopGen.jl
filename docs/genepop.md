## Format

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is tab or space delimited

### Formatting Examples

```tab="loci stacked vertically"
Wasp populations in New York
Locus1
Locus2
Locus3
POP
Oneida_01,  250230 564568 110100
Oneida_02,  252238 568558 100120
Oneida_03,  254230 564558 090100
POP
Newcomb_01,  254230 564558 080100
Newcomb_02,  000230 564558 090080
Newcomb_03,  254230 000000 090100
Newcomb_04,  254230 564000 090120
```

```tab="loci stacked horizontally"
Wasp populations in New York
Locus1,Locus2,Locus3
POP
Oneida_01,  250230 564568 110100
Oneida_02,  252238 568558 100120
Oneida_03,  254230 564558 090100
POP
Newcomb_01,  254230 564558 080100
Newcomb_02,  000230 564558 090080
Newcomb_03,  254230 000000 090100
Newcomb_04,  254230 564000 090120
```

## Import a genepop file as a `PopObj`

!!! warning "Windows users"
    make sure to change your backslashes "\" to forward slashes "/" 

```julia
genepop(infile; digits = 3, popsep = "POP", numpops)

# Example
julia> b = genepop("/data/wasp_hive.gen", digits = 3, popsep = "POP", numpops = 2)
```

### Arguments

- `#!julila infile::String` : path to genepop file, in quotes

### Keyword Arguments

- `#!julila digits::Int64` : the number of digits used to denote an allele (default = 3)
- `#!julila popsep::String` : word that separates populations in `infile` (default: "POP")
- `#!julila numpops::Int64` : number of populations present in `infile` (used for early error checking)

!!! info "Default population names"
    By default, the file reader will assign numbers as population ID's in order of appearance in the genepop file. Use the `popid!` function to rename these with your own population ID's.