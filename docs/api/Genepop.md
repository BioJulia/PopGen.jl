---
id: genepop
title: Genepop.jl
sidebar_label: Genepop.jl
---

### `genepop`
```julia
genepop(infile::String; kwargs...)
```
Load a Genepop format file into memory as a PopData object.
- `infile::String` : path to Genepop file

**Keyword Arguments**
- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent::Bool`   : whether to print file information during import (default: `true`)
- `allow_monomorphic` : whether to keep monomorphic loci in the dataset (default: `false`)

**File must follow standard Genepop formatting**
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default `POP`) must delimit populations
- Sample name is immediately followed by a *comma*
- File is *tab or space delimted* (but not both!)

**Genepop file example:**
```
wasp_hive.gen: Wasp populations in New York \n
Locus1
Locus2
Locus3
pop
Oneida_01,  250230  564568  110100
Oneida_02,  252238  568558  100120
Oneida_03,  254230  564558  090100
pop
Newcomb_01, 254230  564558  080100
Newcomb_02, 000230  564558  090080
Newcomb_03, 254230  000000  090100
Newcomb_04, 254230  564000  090120
```

**Example**
```julia
waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "pop")
```

----

### `popdata2genepop`
```julia
popdata2genepop(data::PopData; filename::String = "output.gen", digits::Int = 3, format::String = "vertical")
```
Writes a `PopData` object to a Genepop-formatted file
- `data`: the `PopData` object you wish to convert to a Genepop file

**Keyword arguments**
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele as (e.g. `(1, 2)` => `001002` for `digits = 3`)
- `format` : a `String` indicating whether loci should be formatted 
    - vertically (`"v"` or `"vertical"`)
    - hortizontally (`"h"`, or `"horizontal"`)
    - Genepop Isolation-By-Distance (`"ibd"`) where each sample is a population with long/lat data prepended

**Example**
```julia
cats = @nancycats;
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
julia> popdata2genepop(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
