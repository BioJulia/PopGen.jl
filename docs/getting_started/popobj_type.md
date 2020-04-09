# PopObj and PopData types
For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called  `PopData`. The struct is defined as:

```julia
struct PopData
	meta::IndexedTable
	loci::IndexedTable
end
```

As you can see, a `PopData` is made up of two IndexedTables (from [JuliaDB.jl](https://github.com/JuliaComputing/JuliaDB.jl)), one called `meta` for sample information (metadata), and the other called `loci` which includes genotype information. This structure allows for easy and convenient access to the fields using dot `.` accessors.. The `meta` and `loci` tables are both specific in their structure, so here is an illustration to help you visualize a `PopData` object:

![PopData](/PopGen.jl/images/PopData.svg)


`PopData` and other custom types introduced in PopGen.jl fall under an AbstractType we call `PopObj`, which is short for "PopGen Object".

::: details pronouncing "PopObj"
It's not super obvious, but we decided to pronounce PopObj as "pop ob" with a silent j because it sounds better than saying "pop obj", but writing it as PopOb looks weird. It's a silly little detail that Pavel seems to care a lot about.
:::

::: danger avoid manual creation!
While it may seem simple enough to create two IndexedTables and make a `PopData` out of them, the structure of `meta` and `loci` are specific, so small mistakes in creating them can create many errors and prevent PopGen from working correctly on your data. Please use the included `csv`, `genepop`, and `vcf` file importers instead.
:::

## Metadata

The `meta` table has 5 specific categories/columns: name, population, ploidy, longitude, latitude. These can be directly accessed with `PopData.meta.columns.colname` where `PopData` is the name of your PopData object, and `colname` is one of the five column names below.

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

## Genotype Information

The genotype information is stored in a separate table called `loci`. This table is rather special in that it is stored in "tidy" format, i.e. one record per row. Storing data this way makes it a lot easier to interrogate the data and write new functions, along with leveraging [JuliaDBMeta.jl](https://github.com/piever/JuliaDBMeta.jl). It also means the table will have as many rows as loci x samples, which can become a lot. To reduce redundant objects inflating object size, the columns name, population, and locus are `CategoricalStrings`  from [CategoricalArrays.jl](https://github.com/JuliaData/CategoricalArrays.jl), which is a memory-saving data structure for long repetitive categorical data. Without using this format, `gulfsharks`, whose source file is 3.2mb, would occupy about 27mb in your RAM! The classes of `.loci` can be directly accessed with `PopData.loci.columns.colname` where `PopData` is the name of your PopData object, and `colname` is one of the four column names below.

### name

`::Vector{CategoricalString}`

The sample name, stored as a `CategoricalString`. Fundamentally, this acts like the `name` column of the `meta` table, except when deleting entries and a few uncommon edge cases.

### population

`::Vector{CategoricalString}`

The population ID associated with that sample, stored as a `CategoricalString`. Fundamentally, this acts like the `population` column of the `meta` table, except when deleting entries and a few uncommon edge cases.

### locus

`::Vector{CategoricalString}`

The locus associated with the genotype, stored as a `CategoricalString`.

### genotype

`::Vector{Union{Missing,Genotype}}`

The genotypes of the `loci` are an array of type `Genotype`, which is an alias for a built-in Julia Tuple type with each value corresponding to an allele (read below to disentangle what that type actually is). For the most part, it looks like this:

```julia tab="genotype example"
[(0,1), (0,0), missing, (1,2)]
```

::: warning immutable genotypes
We use the Tuple type for genotypes of individuals because they are **immutable** (cannot be changed). By the time you're using PopGen.jl, your data should already be filtered and screened. Hand-editing of genotype alleles is **strongly** discouraged, so we outlawed it.
:::

## Viewing PopData

Given the volume of information that can be present in a `PopData`, we recommend `summary()` to summarize/overview the data rather than regurgitate everything on the screen. 

```
julia> a = gulfsharks() ;

julia> summary(a)
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 212
  Number of loci: 2213
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```

## location data

Location data is optional for a `PopData`. There are functions that use location information (e.g. `locations`), but most don't, so it's not a dealbreaker. At present, there are no analyses that utilize location information. 