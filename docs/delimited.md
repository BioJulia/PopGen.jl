## Formatting

- Loci names must be first row
- Individuals names must be first value in row
- Population ID's must be second value in row
- **Optional:** longitude (x) values third value in row, latitude (y) fourth

### Formatting examples

```bash tab="without locations"
Locus1,Locus2,Locus3
sierra_01,1,001001,002002,001001
sierra_02,1,001001,001001,001002
snbarb_03,2,001001,001001,001002
snbarb_02,2,001001,001001,001001
snbarb_03,2,001002,001001,001001
```

```bash tab="with locations"
Locus1,Locus2,Locus3
sierra_01,1,14.1,15.2,001001,002002,001001
sierra_02,1,34.1,26.1,001001,001001,001002
snbarb_03,2,45.1,-11.2,001001,001001,001002
snbarb_02,2,-11.5,11.6,001001,001001,001001
snbarb_03,2,-3.1,43.2,001002,001001,001001
```



## Importing a CSV into Julia as a `PopObj`

```julia
csv(infile; delim, ploidy = 2, location = false)

# Example
a = csv("/data/cali_poppy.csv", delim = ",", ploidy = 2)
```

### Arguments

- `#!julia infile::String` : path to the input file, in quotes

### Keyword Arguments

- `#!julia delim::Union{Char,String,Regex}` : delimiter of the file, can be a string, character, or regex
  - comma: `delim = ","`
  - space: `delim = " "`
  - tab: `delim = "\t"`
  - etc.
- `#!julia ploidy::Int64` : single integer of the ploidy of the samples in the file (default = 2)
  - haploid: `ploidy = 1`
  - diploid: `ploidy = 2`
  - triploid: `ploidy = 3`
  - etc.
- `#!julia location::Bool = false` : true/false of whether location data is present in the file (default = false)