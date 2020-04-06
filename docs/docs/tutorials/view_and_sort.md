# Viewing and Sorting

The functions here help you inspect your `PopData` and pull information from it easily.

## Individuals / Samples

### view individuals' names

```julia
samples(data::PopData)
```

View individual/sample names in a `PopData`. 
:::: tabs card true
::: tab samples
``` julia
julia> samples(sharks)
```
:::
::: tab output
```
212-element Array{String,1}:
 "cc_001" 
 "cc_002" 
 "cc_003" 
 "cc_005" 
 "cc_007" 
 "cc_008" 
 "cc_009" 
 "cc_010" 
 "cc_012" 
 "cc_013" 
 ⋮        
 "seg_023"
 "seg_024"
 "seg_025"
 "seg_026"
 "seg_027"
 "seg_028"
 "seg_029"
 "seg_030"
 "seg_031"
```
:::
::::


## Sorting

```julia
reindex(data::PopData, col::Union{String, Symbol})
```

By default, the `loci` table of `PopData` is sorted by the `name` column. For performance or convenience reasons, you can sort it using any column you want. This will sort the `loci` table of a `PopData` object by column `col` and return new `PopData` object, keeping the original intact. The column names can be Strings or Symbols.

:::: tabs card true
::: tab sort
```julia
sorted_sharks = reindex(sharks, :population)
```
:::
::::


## Display Specific Loci and/or Samples

### Get loci names

```julia
loci(data::PopData)
```

Returns a vector of strings of the loci names in a `PopData`
:::: tabs card true
::: tab loci
```julia
julia> loci(sharks)
```
:::
::: tab output
```
2213-element Array{String,1}:
 "contig_35208"
 "contig_23109"
 "contig_4493" 
 "contig_10742"
 "contig_14898"
 "contig_8483" 
 "contig_8065" 
 "contig_14708"
 "contig_2307" 
 ⋮             
 "contig_24711"
 "contig_18959"
 "contig_43517"
 "contig_27356"
 "contig_475"  
 "contig_19384"
 "contig_22368"
 "contig_2784" 
```
:::
::::

### View genotypes at a locus

```julia
locus(data::PopData, locus::String)
```

Default shows all genotypes for all individuals. Returns a Vector.
:::: tabs card true
::: tab all loci
``` julia
julia> locus(sharks, "contig_2784")
```
:::
::: tab output
```
212-element view(::Array{Union{Missing, Tuple{Int8,Int8}},1}, [2213, 4426, 6639, 8852, 11065, 13278, 15491, 17704, 19917, 22130  …  449239, 451452, 453665, 455878, 458091, 460304, 462517, 464730, 466943, 469156]) with eltype Union{Missing, Tuple{Int8,Int8}}:
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 ⋮      
 (1, 1) 
 (1, 1) 
 (1, 1) 
 (1, 1) 
 missing
 (1, 1) 
 (1, 1) 
 (1, 1) 
```

::: tip returning a view
The `locus` function gives a `view` into the genotype section of the `loci` table.. The text above the  output is isn't anything to be worried about-- it's indicating you are looking at a `view` of a column of the table and that it's not returning a new vector.
:::
::::

## View genotypes by sample (or anything)

This can be done fairly easily using JuliaDBMeta macro `@where`
:::: tabs card true
::: tab single sample
```julia
julia> @where sharks.loci :name == "cc_001"
Table with 2213 rows, 4 columns:
name      population        locus           genotype
────────────────────────────────────────────────────
"cc_001"  "Cape Canaveral"  "contig_35208"  (1, 2)
"cc_001"  "Cape Canaveral"  "contig_23109"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_4493"   (1, 2)
"cc_001"  "Cape Canaveral"  "contig_10742"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_14898"  (1, 2)
"cc_001"  "Cape Canaveral"  "contig_8483"   (1, 1)
"cc_001"  "Cape Canaveral"  "contig_8065"   (1, 1)
"cc_001"  "Cape Canaveral"  "contig_14708"  (1, 1)
⋮
"cc_001"  "Cape Canaveral"  "contig_43517"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_27356"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_475"    (1, 2)
"cc_001"  "Cape Canaveral"  "contig_19384"  (2, 2)
"cc_001"  "Cape Canaveral"  "contig_22368"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_2784"   (1, 1)
```
:::
::: tab multiple samples
```julia
julia> @where sharks.loci :name in ["cc_001", "cc_002"]
Table with 4426 rows, 4 columns:
name      population        locus           genotype
────────────────────────────────────────────────────
"cc_001"  "Cape Canaveral"  "contig_35208"  (1, 2)
"cc_001"  "Cape Canaveral"  "contig_23109"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_4493"   (1, 2)
"cc_001"  "Cape Canaveral"  "contig_10742"  (1, 1)
"cc_001"  "Cape Canaveral"  "contig_14898"  (1, 2)
"cc_001"  "Cape Canaveral"  "contig_8483"   (1, 1)
"cc_001"  "Cape Canaveral"  "contig_8065"   (1, 1)
"cc_001"  "Cape Canaveral"  "contig_14708"  (1, 1)
⋮
"cc_002"  "Cape Canaveral"  "contig_43517"  (1, 1)
"cc_002"  "Cape Canaveral"  "contig_27356"  (1, 1)
"cc_002"  "Cape Canaveral"  "contig_475"    (1, 2)
"cc_002"  "Cape Canaveral"  "contig_19384"  (2, 2)
"cc_002"  "Cape Canaveral"  "contig_22368"  (1, 1)
"cc_002"  "Cape Canaveral"  "contig_2784"   (1, 1)
```
:::
::: tab name and locus
It also means that you can combine different queries with "and" `&&` and "or" `||`. Here is an example of an approach combining a name and locus criteria:

```julia
julia> @where sharks.loci :name in ["cc_001", "cc_002"] && :locus == "contig_2784"
Table with 2 rows, 4 columns:
name      population        locus          genotype
───────────────────────────────────────────────────
"cc_001"  "Cape Canaveral"  "contig_2784"  (1, 1)
"cc_002"  "Cape Canaveral"  "contig_2784"  (1, 1)
```
:::
::::


## Missing Data

```julia
missing(data::PopData; mode::String = "sample")
```

Get missing genotype information in a `PopData`. Specify a `mode` of operation to return a DataFrame corresponding with that missing information type.

| mode     | alternative name | what it does                                                 |
| -------- | ---------------- | ------------------------------------------------------------ |
| `"sample"` | `"individual"`     | returns a count and list of missing loci per individual (default) |
| `"pop"`    | `"population"`     | returns a count of missing genotypes per population          |
| `"locus"`  | `"loci"`           | returns a count of missing genotypes per locus               |
| `"full"`   | `"detailed"`       | returns a count of missing genotypes per locus per population |

:::: tabs card true
::: tab sample
```
julia> missing(sharks)
Table with 212 rows, 2 columns:
name       missing
──────────────────
"cc_001"   124
"cc_002"   94
"cc_003"   100
"cc_005"   0
"cc_007"   2
"cc_008"   1
"cc_009"   2
"cc_010"   1
⋮
"seg_026"  0
"seg_027"  2
"seg_028"  25
"seg_029"  0
"seg_030"  1
"seg_031"  1
```
:::
::: tab pop
```
julia> missing(sharks, mode = "pop")
Table with 7 rows, 2 columns:
population        missing
─────────────────────────
"Florida Keys"    782
"Cape Canaveral"  666
"Mideast Gulf"    379
"Georgia"         744
"Northeast Gulf"  93
"Southeast Gulf"  1504
"South Carolina"  480
```
:::
::: tab locus
```
julia> missing(sharks, mode = "locus")
Table with 2213 rows, 2 columns:
locus           missing
───────────────────────
"contig_35208"  0
"contig_23109"  6
"contig_4493"   3
"contig_10742"  2
"contig_14898"  0
"contig_8483"   0
"contig_8065"   0
"contig_14708"  1
⋮
"contig_43517"  6
"contig_27356"  2
"contig_475"    0
"contig_19384"  5
"contig_22368"  3
"contig_2784"   7
```
:::
::: tab full
```tab
julia> missing(sharks, mode = "full")
Table with 15491 rows, 3 columns:
locus           population        missing
─────────────────────────────────────────
"contig_35208"  "Florida Keys"    0
"contig_35208"  "Cape Canaveral"  0
"contig_35208"  "Mideast Gulf"    0
"contig_35208"  "Georgia"         0
"contig_35208"  "Northeast Gulf"  0
"contig_35208"  "Southeast Gulf"  0
"contig_35208"  "South Carolina"  0
"contig_23109"  "Florida Keys"    0
⋮
"contig_2784"   "Cape Canaveral"  0
"contig_2784"   "Mideast Gulf"    1
"contig_2784"   "Georgia"         1
"contig_2784"   "Northeast Gulf"  1
"contig_2784"   "Southeast Gulf"  1
"contig_2784"   "South Carolina"  0
```
:::
::::
::: tip alternative names
Each `mode` has an extra synonymous (alternative) name just because we can and want you to have the option of more explicitly legible code. If you get the `mode` wrong, it will let you know with an error message and run the default `"sample"` mode anyway.
:::