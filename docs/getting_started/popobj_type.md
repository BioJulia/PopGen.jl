---
id: types
title: PopObj and PopData types
sidebar_label: PopObj and PopData types
---

For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called  `PopData`. The struct is defined as:

```julia
struct PopData
	meta::DataFrame
	loci::DataFrame
end
```

As you can see, a `PopData` is made up of two DataFrames, one called `meta` for sample information (metadata), and the other called `loci` which includes genotype information. This structure allows for easy and convenient access to the fields using dot `.` accessors.. The `meta` and `loci` tables are both specific in their structure, so here is an illustration to help you visualize a `PopData` object:

![PopData](/img/PopData.svg)

`PopData` and other custom types introduced in PopGen.jl fall under an AbstractType we call `PopObj`, which is short for "PopGen Object".

:::note pronouncing "PopObj"
It's not super obvious, but we decided to pronounce PopObj as "pop ob" with a silent j because it sounds better than saying "pop obj", but writing it as PopOb looks weird. It's a silly little detail that Pavel seems to care a lot about.
:::

:::danger avoid manual creation!
While it may seem simple enough to create two DataFrames and make a `PopData` out of them, the structure of `meta` and `loci` are specific, so small mistakes in creating them can create many errors and prevent PopGen from working correctly on your data. Please use the included `csv`, `genepop`, and `vcf` file importers instead.
:::

## Metadata

The `meta` table has 5 specific categories/columns: name, population, ploidy, longitude, latitude. These can be directly accessed with `PopData.meta.colname` where `PopData` is the name of your PopData object, and `colname` is one of the five column names below.

### name

`::Vector{String}`

The individual/sample names

```julia
["ind_001", "ind_002", "ind_003"]
```

### population

`::Vector{String}`

The individual/sample population ID's

```julia
["borneo", "borneo", "new jersey"]
```

### ploidy

`::Vector{Int8}`

The ploidy of the samples

```julia
[2, 2, 2]
```

### longitude

`::Vector{Union{Missing,Float32}}`

latitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4.42]
```

### latitude

`::Vector{Union{Missing,Float64}}`

longitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```


:::note location data
Location data is optional for `PopData`. There are functions that use location information (e.g. `locations`), but most don't, so it's not a dealbreaker. At present, there are no analyses that utilize location information. 
:::


## Genotype Information

The genotype information is stored in a separate table called `loci`. This table is rather special in that it is stored in "tidy" format, i.e. one record per row. Storing data this way makes it a lot easier to interrogate the data and write new functions. It also means the table will have as many rows as loci x samples, which can become a lot. To reduce redundant objects inflating object size, the columns name, population, and locus are each a special type of compressed vector from [PooledArrays.jl](https://github.com/JuliaData/PooledArrays.jl), which is a memory-saving data structure for long repetitive categorical data. Without using this format, `gulfsharks`, whose source file is 3.2mb, would occupy about 27mb in your RAM! The classes of `.loci` can be directly accessed with `PopData.loci.colname` where `PopData` is the name of your PopData object, and `colname` is one of the four column names below. For clarity, the columns will be represented below as though they are regular vectors.

### name

`::Vector{String}`

The sample name, stored as a `PooledArray` of eltype `String`. Fundamentally, this acts like the `name` column of the `meta` table, except when deleting entries and a few uncommon edge cases.

### population

`::Vector{String}`

The population ID associated with that sample,stored as a `PooledArray` of eltype `String`. Fundamentally, this acts like the `population` column of the `meta` table, except when deleting entries and a few uncommon edge cases.

### locus

`::Vector{String}`

The locus associated with the genotype, stored as a `PooledArray` of eltype `String`.

### genotype

`::Vector{Union{Missing,Genotype}}`

The genotypes of the `loci` are an array of type `Genotype`, which is [an alias](/getting_started/other_types.md) for a built-in Julia Tuple type with each value corresponding to an allele. For the most part, it looks like this:

```julia
[(0,1), (0,0), missing, (1,2)]
```

:::caution immutable genotypes
We use the Tuple type for genotypes of individuals because they are **immutable** (cannot be changed). By the time you're using PopGen.jl, your data should already be filtered and screened. Hand-editing of genotype alleles is **strongly** discouraged, so we outlawed it.
:::

------

## Acknowledgements
A *lot* of what's possible in PopGen.jl is thanks to the tireless work of the contributors and maintainers of [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl). It's no small task to come up with and maintain a robust, performant, and sensible tabular data type, and they deserve so much credit for it. 