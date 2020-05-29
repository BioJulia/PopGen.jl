---
id: viewsort
title: Viewing and sorting
sidebar_label: Viewing and sorting
---

The functions here help you inspect your `PopData` and pull information from it easily.

## Individuals / Samples

### view individuals' names

```julia
samples(data::PopData)
```

View individual/sample names in a `PopData`. 
:::: tabs card stretch
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


## Display Specific Loci and/or Samples

### Get loci names

```julia
loci(data::PopData)
```

Returns a vector of strings of the loci names in a `PopData`
:::: tabs card stretch
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
:::: tabs card stretch
::: tab all loci

``` julia
julia> locus(sharks, "contig_2784")
```

:::
::: tab output

```
212-element Array{Union{Missing, Tuple{Int8,Int8}},1}:
 (1, 1)
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

::::

## View genotypes by sample (or anything)

This can be done fairly easily using DataFramesMeta macro `@where`
:::: tabs card stretch
::: tab single sample

```julia
julia> @where(sharks.loci, :name .== "cc_001")
2213×4 DataFrame
│ Row  │ name   │ population     │ locus        │ genotype │
│      │ Cat…   │ Categorical…   │ Categorical… │ Tuple…?  │
├──────┼────────┼────────────────┼──────────────┼──────────┤
│ 1    │ cc_001 │ Cape Canaveral │ contig_35208 │ (1, 2)   │
│ 2    │ cc_001 │ Cape Canaveral │ contig_23109 │ (1, 1)   │
│ 3    │ cc_001 │ Cape Canaveral │ contig_4493  │ (1, 2)   │
│ 4    │ cc_001 │ Cape Canaveral │ contig_10742 │ (1, 1)   │
│ 5    │ cc_001 │ Cape Canaveral │ contig_14898 │ (1, 2)   │
│ 6    │ cc_001 │ Cape Canaveral │ contig_8483  │ (1, 1)   │
⋮
│ 2207 │ cc_001 │ Cape Canaveral │ contig_18959 │ (1, 2)   │
│ 2208 │ cc_001 │ Cape Canaveral │ contig_43517 │ (1, 1)   │
│ 2209 │ cc_001 │ Cape Canaveral │ contig_27356 │ (1, 1)   │
│ 2210 │ cc_001 │ Cape Canaveral │ contig_475   │ (1, 2)   │
│ 2211 │ cc_001 │ Cape Canaveral │ contig_19384 │ (2, 2)   │
│ 2212 │ cc_001 │ Cape Canaveral │ contig_22368 │ (1, 1)   │
│ 2213 │ cc_001 │ Cape Canaveral │ contig_2784  │ (1, 1)   │
```

:::
::: tab multiple samples

```julia
julia> @where(sharks.loci, :name .∈ Ref(["cc_001", "cc_002"]))
4426×4 DataFrame
│ Row  │ name   │ population     │ locus        │ genotype │
│      │ Cat…   │ Categorical…   │ Categorical… │ Tuple…?  │
├──────┼────────┼────────────────┼──────────────┼──────────┤
│ 1    │ cc_001 │ Cape Canaveral │ contig_35208 │ (1, 2)   │
│ 2    │ cc_002 │ Cape Canaveral │ contig_35208 │ (1, 2)   │
│ 3    │ cc_001 │ Cape Canaveral │ contig_23109 │ (1, 1)   │
│ 4    │ cc_002 │ Cape Canaveral │ contig_23109 │ (1, 2)   │
│ 5    │ cc_001 │ Cape Canaveral │ contig_4493  │ (1, 2)   │
│ 6    │ cc_002 │ Cape Canaveral │ contig_4493  │ (1, 1)   │
⋮
│ 4420 │ cc_002 │ Cape Canaveral │ contig_475   │ (1, 2)   │
│ 4421 │ cc_001 │ Cape Canaveral │ contig_19384 │ (2, 2)   │
│ 4422 │ cc_002 │ Cape Canaveral │ contig_19384 │ (2, 2)   │
│ 4423 │ cc_001 │ Cape Canaveral │ contig_22368 │ (1, 1)   │
│ 4424 │ cc_002 │ Cape Canaveral │ contig_22368 │ (1, 1)   │
│ 4425 │ cc_001 │ Cape Canaveral │ contig_2784  │ (1, 1)   │
│ 4426 │ cc_002 │ Cape Canaveral │ contig_2784  │ (1, 1)   │
```

:::
::: tab name and locus
It also means that you can combine different queries with "and" `&&` and "or" `||`. Here is an example of an approach combining a name and locus criteria:

```julia
julia> @where(sharks.loci, :name .∈ Ref(["cc_001", "cc_002"]), :locus .== "contig_2784")
2×4 DataFrame
│ Row │ name   │ population     │ locus       │ genotype │
│     │ Cat…   │ Categorical…   │ Cat…        │ Tuple…?  │
├─────┼────────┼────────────────┼─────────────┼──────────┤
│ 1   │ cc_001 │ Cape Canaveral │ contig_2784 │ (1, 1)   │
│ 2   │ cc_002 │ Cape Canaveral │ contig_2784 │ (1, 1)   │
```

:::
::::


## Missing Data

```julia
missing(data::PopData; by::String = "sample")
```

Get missing genotype information in a `PopData`. Specify a mode of operation using the `by =` keyword to return a table corresponding with that missing information type.

|     by     | alternative name | what it does                                                 |
| :--------: | :--------------: | ------------------------------------------------------------ |
| `"sample"` |  `"individual"`  | returns a count of missing loci per individual (default)     |
|  `"pop"`   |  `"population"`  | returns a count of missing genotypes per population          |
| `"locus"`  |     `"loci"`     | returns a count of missing genotypes per locus               |
|  `"full"`  |   `"detailed"`   | returns a count of missing genotypes per locus per population |

:::: tabs card stretch
::: tab sample

```
julia> missing(sharks)
212×2 DataFrame
│ Row │ name    │ missing │
│     │ Cat…    │ Int64   │
├─────┼─────────┼─────────┤
│ 1   │ cc_001  │ 124     │
│ 2   │ cc_002  │ 94      │
│ 3   │ cc_003  │ 100     │
│ 4   │ cc_005  │ 0       │
│ 5   │ cc_007  │ 2       │
│ 6   │ cc_008  │ 1       │
⋮
│ 206 │ seg_025 │ 0       │
│ 207 │ seg_026 │ 0       │
│ 208 │ seg_027 │ 2       │
│ 209 │ seg_028 │ 25      │
│ 210 │ seg_029 │ 0       │
│ 211 │ seg_030 │ 1       │
│ 212 │ seg_031 │ 1       │
```

:::
::: tab pop

```
julia> missing(sharks, by = "pop")
7×2 DataFrame
│ Row │ population     │ missing │
│     │ Categorical…   │ Int64   │
├─────┼────────────────┼─────────┤
│ 1   │ Florida Keys   │ 1246    │
│ 2   │ Cape Canaveral │ 666     │
│ 3   │ Mideast Gulf   │ 99      │
│ 4   │ Georgia        │ 425     │
│ 5   │ Northeast Gulf │ 474     │
│ 6   │ Southeast Gulf │ 1504    │
│ 7   │ South Carolina │ 234     │
```

:::
::: tab locus

```
julia> missing(sharks, by = "locus")
2213×2 DataFrame
│ Row  │ locus        │ missing │
│      │ Categorical… │ Int64   │
├──────┼──────────────┼─────────┤
│ 1    │ contig_35208 │ 0       │
│ 2    │ contig_23109 │ 6       │
│ 3    │ contig_4493  │ 3       │
│ 4    │ contig_10742 │ 2       │
│ 5    │ contig_14898 │ 0       │
│ 6    │ contig_8483  │ 0       │
⋮
│ 2207 │ contig_18959 │ 0       │
│ 2208 │ contig_43517 │ 6       │
│ 2209 │ contig_27356 │ 2       │
│ 2210 │ contig_475   │ 0       │
│ 2211 │ contig_19384 │ 5       │
│ 2212 │ contig_22368 │ 3       │
│ 2213 │ contig_2784  │ 7       │
```

:::
::: tab full

```
julia> missing(sharks, by = "full")
15491×3 DataFrame
│ Row   │ locus        │ population     │ missing │
│       │ Categorical… │ Categorical…   │ Int64   │
├───────┼──────────────┼────────────────┼─────────┤
│ 1     │ contig_35208 │ Cape Canaveral │ 0       │
│ 2     │ contig_35208 │ Georgia        │ 0       │
│ 3     │ contig_35208 │ South Carolina │ 0       │
│ 4     │ contig_35208 │ Florida Keys   │ 0       │
│ 5     │ contig_35208 │ Mideast Gulf   │ 0       │
│ 6     │ contig_35208 │ Northeast Gulf │ 0       │
⋮
│ 15485 │ contig_2784  │ Cape Canaveral │ 0       │
│ 15486 │ contig_2784  │ Georgia        │ 2       │
│ 15487 │ contig_2784  │ South Carolina │ 1       │
│ 15488 │ contig_2784  │ Florida Keys   │ 2       │
│ 15489 │ contig_2784  │ Mideast Gulf   │ 1       │
│ 15490 │ contig_2784  │ Northeast Gulf │ 0       │
│ 15491 │ contig_2784  │ Southeast Gulf │ 1       │
```

:::
::::
:::tip alternative names
Each mode of operation has an extra synonymous (alternative) name just because we can and want you to have the option of more explicitly legible code. If you get the `by = `  argument wrong, it will let you know with an error message and run the default `"sample"` mode anyway.
:::