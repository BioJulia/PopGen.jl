For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called a `PopObj`. The struct is defined as:

```julia
mutable struct PopObj
	samples::DataFrame
	loci::DataFrame
end
```

!!! info "pronouncing "PopObj" "
    If you haven't already guessed, `PopObj` is a combination of the words PopGen and Object. PopObj is pronounced "pop ob" with a silent j because it rolls of the tongue better, but writing it as PopOb looks weird. 
    
    Yes, I have lost sleep over this detail.    
    - Pavel

As you can see, a `PopObj` is made up of two dataframes, one for sample information, the other for genotype information. This structure allows for easy and convenient access to the fields using dot `.` accessors.

!!! failure "avoid manual creation"
    While it may seem simple enough to create two dataframes and make a `PopObj` out of them, the structure of `samples` and `loci` are specific, so small mistakes in creating them can create many errors and prevent PopGen from working correctly on your data. Please use the included `csv`, `genepop`, and `vcf` file importers instead.

## samples

The `samples` dataframe has 5 specific categories: name, population, ploidy, latitude, longitude.

### `samples.name` 

`::Array{String,1}`

The individual/sample names

```julia
["ind_001", "ind_002", "ind_003"]
```

### `samples.population`

`::Array{String,1}`

The individual/sample population ID's

```julia
["borneo", "borneo", "new jersey"]
```

### `samples.ploidy`

`::Array{Int8,1}`

The ploidy of the samples

```julia
[2, 2, 2]
```

### `samples.latitude`

`::Array{Union{Int64,Float64},1}`

latitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```

### `samples.longitude`

`::Array{Union{Int64,Float64},1}`

longitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```

## loci

The genotype information is stored in a separate dataframe called `loci`, where each column is named for the locus it represents. This makes for easy and obvious accessing by calling `PopObj.loci.locusname`.  To view the loci names, use the convenient `loci_names` command.

### genotypes 

`::Array{Tuple{Int16,...},1}`

The genotypes of the `loci` are an array of tuples, with each value corresponding to an allele. The length of the tuple will vary based on the ploidy of the sample, therefor the `type` shown above is conceptually accurate, but computationally incorrect.

```
[(0,1),(0,0),(1,2)]
```

!!! important
    We use the tuple type for genotypes of individuals because they are **immutable** (cannot be changed). By the time you're using `PopGen.jl`, your data should already be filtered and screened. Hand-editing of genotype values is **strongly** discouraged, so we outlawed it outright.

## Viewing a PopObj

Given the volume of information that can be present in a `PopObj`, we recommend `summary()` to summarize/overview the data rather than regurgitate everything on the screen. 

```
julia> a = gulfsharks() ;

julia> summary(a)
Object of type PopObj:
No location data provided

Number of individuals: 212
["cca_001", "cca_002", "cca_003"] … ["seg_029", "seg_030", "seg_031"]

Number of loci: 2213
["Contig_35208", "Contig_23109", "Contig_4493"] … ["Contig_19384", "Contig_22368", "Contig_2784"]

Ploidy: 2
Number of populations: 7

Population names and counts:
7×2 DataFrames.DataFrame
│ Row │ population       │ count │
│     │ Categorical…⍰    │ Int32 │
├─────┼──────────────────┼───────┤
│ 1   │ "Cape Canaveral" │ 21    │
│ 2   │ "Georgia"        │ 30    │
│ 3   │ "South Carolina" │ 28    │
│ 4   │ "Florida Keys"   │ 65    │
│ 5   │ "Mideast Gulf"   │ 28    │
│ 6   │ "Northeast Gulf" │ 20    │
│ 7   │ "Southeast Gulf" │ 20    │

Available .samples fields: .name, .population, .ploidy, .longitude, .latitude
```



!!! info "the secret "PopOpt" type"
    There is a complementary type to the `PopObj` called the `PopOpt` ("PopGen Optimized"), which isn't exported, but callable with `PopGen.PopOpt(x::PopObj)`. This is a behind-the-scenes immutable version of a `PopObj` that exists to boost performance and efficiency. The first thing some PopGen.jl commands do is make a temporary `PopOpt` copy of you `PopObj` and use that for indexing, sorting, etc. Using this method allows us to speed up runtime 2x-10x, substantially reduce RAM usage for commands, and still give you the flexibility to augment your `PopObj` as needed. You don't need to know this bit of trivia to use PopGen.jl, but it may be useful if you plan on writing your own functions. 