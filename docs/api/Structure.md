---
id: structure
title: Structure.jl
sidebar_label: Structure.jl
---

### `phase_structure`
```julia
phase_structure(datatype::DataType, args...)
```
Takes a DataType (such as `Int8`) and a series of integers to return
a sorted Tuple of those integers converted to that DataType. i.e. takes
a series of alleles and returns a genotype. Returns `missing` if args are
`missing`. Used internally in PopGen.structure file reader.

#### Example
```
phase_structure(Int8, 1,2,3,4,3,4,6,1)
(1, 1, 2, 3, 3, 4, 4, 6)

phase_structure(Int16, missing, missing)
missing
```
----
### `structure`
```julia
    structure(infile::String; kwargs...)
```
Load a Structure format file into memory as a PopData object.

- `infile::String` : path to Structure file

#### Keyword Arguments
- `extracols::Integer`: how many additional optional columns there are beyond Stucture's POPDATA the reader needs to ignore (default: `0`)
    - these include POPFLAG, LOCDATA, or anything else you might have added
- `extrarows::Integer` : how many additional optional rows there are beyond the first row of locus names (default: `0`)
- `missingval::String`  : the value used to identify missing values in the data (default: `"-9"`)
- `silent::Bool`   : whether to print file information during import (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)
- `faststructure::Bool`: whether the file is fastStructure format (default: `false`)

#### File must follow this Structure format:
- the file is `tab` or `space` delimited **but not both**
- first row is locus names separated by the delimiter
    - leading/trailing whitespaces are tolerated
    - optional rows allowed **after** the locus names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
    - multi-column variant not supported
- first data column is sample name
- second data column is population ID
    - optional columns allowed **after** the population ID (2nd) column
- remaining columns are the genotype for that individual for that locus

#### Structure file example:
```
locus_1	locus_2	locus_3	locus_4	locus_5
walnut_01	1	-9	145	66	0	92
walnut_01	1	-9	-9	64	0	94
walnut_02	1	106	142	68	1	92
walnut_02	1	106	148	64	0	94
walnut_03	2	110	145	-9	0	92
walnut_03	2	110	148	66	1	-9
```

#### fastStructure file format:
- the file is `tab` or `space` delimited **but not both**
- no first row of loci names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
- first data column is sample name
- second data column is population ID
- remaining columns are the genotype for that individual for that locus
- usually, first 6 colums are empty (but not necessary)
- **no** extra rows or columns.
#### fastStructure file example:
```
chestnut_01	1	-9	145	66	0	92
chestnut_01	1	-9	-9	64	0	94
chestnut_02	1	106	142	68	1	92
chestnut_02	1	106	148	64	0	94
chestnut_03	2	110	145	-9	0	92
chestnut_03	2	110	148	66	1	-9
```

#### Example
```
walnuts = structure("juglans_nigra.str", extracols = 0, extrarows = 0)
```
----

### `popdata2structure`
```julia
popdata2structure(data::PopData; filename::String, faststructure::Bool, delim::String)
```
Write a `PopData` object to a Stucture format file
- `data`: the `PopData` object you wish to convert to a Structure file
### keyword arguments
- `filename`: a `String` of the output filename
- `delim` : a `String` of either `"tab"` or `"space"` indicating the delimiter (default: `"tab"`)
- `faststructure`: true/false of whether the output should be formatted for fastStructure (default: `false`)

```
cats = @nancycats;
fewer_cats = omit(cats, name = samples(cats)[1:10]);
popdata2structure(fewer_cats, filename = "filtered_nancycats.str", faststructure = true)
```