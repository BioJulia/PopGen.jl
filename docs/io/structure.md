---
id: structure
title: Structure format
sidebar_label: Structure format
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

:::info
More often than not, your Structure file was created by a conversion from another format. While PopGen.jl offers a Structure file reader, we generally recommend using whatever previous format it was in because the Structure reader has more specific format requirements than the other readers, which can cause unneeded frustration. Additionally, fewer data conversions mean less chance of conversion errors occuring. 
:::

## Import a Structure file as `PopData`


```julia
structure(infile::String; kwargs...)
```

### Arguments
- `infile::String` : path to Structure file

### Keyword Arguments
- `extracols::Integer`: how many additional optional columns there are beyond Stucture's POPDATA the reader needs to ignore (default: `0`)
    - these include POPFLAG, LOCDATA, or anything else you might have added
- `extrarows::Integer` : how many additional optional rows there are beyond the first row of locus names (default: `0`)
- `missingval::String`  : the value used to identify missing values in the data (default: `"-9"`)
- `silent::Bool`   : whether to print file information during import (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)
- `faststructure::Bool`: whether the file is fastStructure format (default: `false`)

### File formatting
Structure files are not an ideal format because there is a bit too much wiggle room in the specifications that are later cleaned up with a config file when running the software. As such, PopGen.jl requires somewhat more specificity in how the files are formatted for things to work correctly:

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'Stucture', value: 's', },
    { label: 'fastStructre', value: 'f', },
  ]
}>
<TabItem value="s">

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
#### Example
```
walnuts = structure("juglans_nigra.str", extracols = 0, extrarows = 0)
```

</TabItem>
<TabItem value="f">

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
chestnuts = structure("castanea_dentata.str", faststructure = true)
```

</TabItem>
</Tabs>


## Writing to a Structure file
All file writing options can be performed using `write_to()`, which calls `popdata2structure` when writing to a Structure/fastStructure file.

```julia
popdata2structure(data::PopData; filename::String, faststructure::Bool, delim::String)
```
Write a `PopData` object to a Stucture format file
### Arguments
- `data`: the `PopData` object you wish to convert to a Structure file

### Keyword Arguments
- `filename`: a `String` of the output filename
- `delim` : a `String` of either `"tab"` or `"space"` indicating the delimiter (default: `"tab"`)
- `faststructure`: true/false of whether the output should be formatted for fastStructure (default: `false`)

#### Example
```
cats = @nancycats;
fewer_cats = omit(cats, name = samples(cats)[1:10]);
popdata2structure(fewer_cats, filename = "filtered_nancycats.str", faststructure = true)
```
