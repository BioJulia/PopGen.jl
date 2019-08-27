# Importing Data

At present, `PopGen.jl` has two file importers: one for CSV (or other delimited) files and another for genepop formatted files. 



## CSV

- Loci names must be first row
- Individuals names must be first value in row
- Population ID's must be second value in row
- **Optional:** longitude (x) values third value in row, latitude (y) fourth

**Formatting examples**

``` bash tab="without locations"
Locus1,Locus2,Locus3
sierra_01,1,001001,002002,001001
sierra_02,1,001001,001001,001002
snbarb_03,2,001001,001001,001002
snbarb_02,2,001001,001001,001001
snbarb_03,2,001002,001001,001001
```

``` bash tab="with locations"
Locus1,Locus2,Locus3
sierra_01,1,14.1,15.2,001001,002002,001001
sierra_02,1,34.1,26.1,001001,001001,001002
snbarb_03,2,45.1,-11.2,001001,001001,001002
snbarb_02,2,-11.5,11.6,001001,001001,001001
snbarb_03,2,-3.1,43.2,001002,001001,001001
```

To import a CSV into Julia as a `PopObj`:

```julia
csv(infile; delim, ploidy = 2, location = false)

# Example
a = csv("/data/cali_poppy.csv", delim = ",", ploidy = 2)
```

Arguments:

-  `infile::String` : path to the input file, in quotes

Keyword Arguments:

- `delim::Union{Char,String,Regex}` : delimiter of the file, can be a string, character, or regex
    - comma: `delim = ","`
    - space: `delim = " "`
    - tab: `delim = "\t"`
    - etc.
- `ploidy::Int64` : single integer of the ploidy of the samples in the file (default = 2)
    - haploid: `ploidy = 1`
    - diploid: `ploidy = 2`
    - triploid: `ploidy = 3`
    - etc.
- `location::Bool = false` : true/false of whether location data is present in the file (default = false)



## Genepop

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is tab or space delimted

**Formatting Examples**

```  tab="loci stacked vertically"
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

```  tab="loci stacked horizontally"
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

To import a genepop file into Julia as a `PopObj`:

```julia
genepop(infile; ploidy = 2, popsep = "POP", numpop)

# Example
b = genepop("/data/wasp_hive.gen", ploidy = 2, popsep = "POP", numpops = 2)
```

Arguments:

- `infile::String` : path to genepop file, in quotes

Keyword Arguments:

- `ploidy::Int64` : single integer of the ploidy of the samples in the file (default = 2)
    - haploid: `ploidy = 1`
    - diploid: `ploidy = 2`
    - triploid: `ploidy = 3`
    - etc.
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `numpops::Int64` : number of populations present in `infile` (used for early error checking)

!!! note
    By default, the file reader will assign numbers as population ID's in order of appearance in the genepop file. Use the `popid!` function to rename these with your own population ID's.

VCF
---------------------

On the To-Do list!