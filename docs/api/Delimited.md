---
id: delimited
title: Delimited.jl
sidebar_label: Delimited.jl
---

### `delimited`
```julia
delimited(infile::String; delim::Union{Char,String,Regex} = "auto", digits::Int64 = 3, silent::Bool = false)
```
Load a delimited-type file into memory as a PopData object. *There should be no empty cells
in your file*
#### Arguments
- `infile` : path to file
#### Keyword Arguments
- `delim` : delimiter characters. By default uses auto-parsing of `CSV.File`
- `digits` : number of digits denoting each allele (default: `3`)
- `diploid`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent`   : whether to print file information during import (default: `true`)
- `allow_monomorphic` : whether to keep monomorphic loci in the dataset (default: `false`)

#### File formatting:
- The first row is column names (names don't matter)
- The columns must be in this order:
    1. sample name
    2. population name
    3. longitude
    4. latitude
    5. locus_1 genotypes
    6. locus_2 genotypes
    7. etc...
#### Missing data
##### Genotypes
Missing genotypes can be formatted as all-zeros `000000`, left empty, or negative-nine `-9`
##### Location data
If location data is missing for a sample (which is ok!), make sure the value is
blank, otherwise there will be transcription errors! (look at line 3 in the example below)
#### Example
```
lizardsCA = delimited("CA_lizards.csv", digits = 3);
```
#### Formatting example

```
name,population,long,lat,Locus1,Locus2,Locus3
sierra_01,mountain,11.11,-22.22,001001,-9,001001
sierra_02,mountain,11.12,-22.21,001001,001001,001002
snbarb_01,coast,,,001001,001001,001002
snbarb_02,coast,11.14,-22.24,001001,001001,001001
snbarb_03,coast,11.15,,001002,001001,001001
```

**Example**
```julia
lizardsCA = delimited("CA_lizards.csv", digits = 3);
```

----

### `popdata2delimtied`
```julia
popdata2delimited(data::PopData; filename::String, delim::String = ",", digits::Integer = 3, format::String = "wide")
```
Write PopData to a text-delimited file. 
#### Keyword Arguments
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele as (e.g. `(1, 2)` => `001002` for `digits = 3`)
- `format` : a `String` indicating whether to output in`"wide"` or `"long"` (aka `"tidy"`) format 
- `delim` : the `String` delimiter to use for writing the file. 

#### Example
```julia
cats = @nancycats;
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
popdata2delimited(fewer_cats, filename = "filtered_nancycats.csv", digits = 3, format = "wide", delim = " ")
```
