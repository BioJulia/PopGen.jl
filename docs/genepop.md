## Format

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is tab or space delimted

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

## Import a genepop file into Julia as a `PopObj`

```julia
genepop(infile; ploidy = 2, popsep = "POP", numpop)

# Example
b = genepop("/data/wasp_hive.gen", ploidy = 2, popsep = "POP", numpops = 2)
```

### Arguments

- `#!julila infile::String` : path to genepop file, in quotes

### Keyword Arguments

- `#!julila ploidy::Int64` : single integer of the ploidy of the samples in the file (default = 2)
  - haploid: `ploidy = 1`
  - diploid: `ploidy = 2`
  - triploid: `ploidy = 3`
  - etc.
- `#!julila popsep::String` : word that separates populations in `infile` (default: "POP")
- `#!julila numpops::Int64` : number of populations present in `infile` (used for early error checking)

!!! note
    By default, the file reader will assign numbers as population ID's in order of appearance in the genepop file. Use the `popid!` function to rename these with your own population ID's.