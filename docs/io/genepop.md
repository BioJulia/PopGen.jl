---
id: genepop
title: Genepop format
sidebar_label: Genepop format
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## Import a genepop file as `PopData`

```julia
genepop(infile; kwargs...)
```

### Arguments

- `infile::String` : path to genepop file, in quotes

### Keyword Arguments

- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent::Bool`   : whether to print file information during import (default: `false`)

:::info population names
By default, the file reader will assign numbers as population ID's (as Strings) in order of appearance in the genepop file. Use the `populations!` function to rename these with your own population ID's.
:::

### Example
```julia
julia> wasp_data = genepop("/data/wasp_hive.gen", digits = 3, popsep = "POP")
```

### Format

Files must follow standard Genepop formatting:

- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular and consistent keyword must delimit populations
- **Must** be the same word each time and not a unique population name
- File is **tab** delimited or **space** delimited, but not both

<Tabs
  block={true}
  defaultValue="v"
  values={[
    { label: 'genepop w/loci stacked vertically', value: 'v', },
    { label: 'genepop w/loci stacked horizontally', value: 'h', },
  ]
}>
<TabItem value="v">

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

</TabItem>
<TabItem value="h">

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

</TabItem>
</Tabs>

## Writing to a Genepop file
All file writing options can be performed using `write_to()`, which calls `popdata2genpop` when writing to a Genepop file.
```julia
popdata2genepop(data::PopData; filename::String = "output.gen", digits::Int = 3, format::String = "vertical", miss::Int = 0)
```
Writes a `PopData` object to a Genepop-formatted file.
### Arguments
- `data`: the `PopData` object you wish to convert to a Genepop file
### Keyword arguments
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele
  -  e.g. `(1, 2)` => `001002` for `digits = 3`
- `format` : a `String` indicating whether loci should be formatted 
  - vertically (`"v"` or `"vertical"`)
  - hortizontally (`"h"`, or `"horizontal"`)
  - Genepop Isolation-By-Distance (`"ibd"`) where each sample is a population with long/lat data prepended
- `miss` : an `Integer` for how you would like missing values written 
  - `0` : As a genotype represented as a number of zeroes equal to `digits × ploidy` like `000000` (default) 
  - `-9` : As a single value `-9`

### Example
```julia
cats = @nancycats;
fewer_cats = omit(cats, names = samples(cats)[1:10]);
julia> popdata2genepop(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```

----

## Acknowledgements

The original implementations of the importing parser were written using only Base Julia, and while the speed was fantastic, the memory footprint involved seemed unusually high (~650mb RAM to parse `gulfsharks`, which is only 3.2mb in size). However, thanks to the efforts of [CSV.jl](https://github.com/JuliaData/CSV.jl), we leverage that package to preserve the speed and reduce the memory footprint quite a bit!
