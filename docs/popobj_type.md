For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called a `PopObj`. The struct is defined as:

```julia
struct PopObj
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

`::Vector{Tuple{Int16,...}}`

The genotypes of the `loci` are an array of tuples, with each value corresponding to an allele. The length of the tuple will vary based on the ploidy of the sample, therefor the `type` shown above is conceptually accurate, but computationally incorrect.

```
[(0,1),(0,0),(1,2)]
```

!!! important
    We use the tuple type for genotypes of individuals because they are **immutable** (cannot be changed). By the time you're using `PopGen.jl`, your data should already be filtered and screened. Hand-editing of genotype values is **strongly** discouraged, so we outlawed it outright.

## viewing a PopObj

Given the volume of information that can be present in a `PopObj`, we recommend `summary()` to summarize/overview the data rather than regurgitate everything on the screen. 

```
julia> a = gulfsharks() ;

julia> summary(a)
 Object of type PopObj
 Marker type: SNP
 Ploidy: 2

 Number of individuals: 212
 Number of loci: 2213
 Longitude: present with 0 missing
 Latitude: present with 0 missing

 Population names and counts:
7×2 DataFrame
│ Row │ population     │ count │
│     │ Union…         │ Int64 │
├─────┼────────────────┼───────┤
│ 1   │ Cape Canaveral │ 21    │
│ 2   │ Georgia        │ 30    │
│ 3   │ South Carolina │ 28    │
│ 4   │ Florida Keys   │ 65    │
│ 5   │ Mideast Gulf   │ 28    │
│ 6   │ Northeast Gulf │ 20    │
│ 7   │ Southeast Gulf │ 20    │
```

## location data

Location data is optional for a `PopObj`. There are functions that use location information (e.g. `locations`and `plot_locations`), but most don't, so it's not a dealbreaker.