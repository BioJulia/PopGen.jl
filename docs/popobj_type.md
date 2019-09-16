For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called a `PopObj`. The struct is defined as:

```julia
mutable struct PopObj
<<<<<<< HEAD
	samples::DataFrame
	loci::DataFrame
=======
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genotypes::Dict
    longitude::Array{Union{Int64,Float64},1}
    latitude::Array{Union{Int64,Float64},1}
>>>>>>> master
end
```

!!! info "pronouncing "PopObj" "
    If you haven't already guessed, `PopObj` is a combination of the words Population and Object. PopObj is pronounced "pop ob" with a silent j because it rolls of the tongue better, but writing it as PopOb looks weird. 
    
    Yes, I have lost sleep over this detail.    
    - Pavel

<<<<<<< HEAD


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
=======
### ind 

**type**: 1-dimension array of strings `::Array{String,1}`

The individual/sample names

```julia
["ind_001", "ind_002", "ind_003"]
```



### popid

**type**: 1-dimension array of strings `::Array{String,1}`

The individual/sample population ID's (in the same order)

```julia
["borneo", "borneo", "new jersey"]
>>>>>>> master
```



<<<<<<< HEAD
### `samples.population`

`::Array{String,1}`

The individual/sample population ID's

```julia
["borneo", "borneo", "new jersey"]
```

### `samples.ploidy`

`::Array{Int8,1}`
=======
### loci

**type**: 1-dimension array of strings `::Array{String,1}`

The name of the loci, as an array of strings

```julia
["locus_001", "locus_2","super-awesome-locus-3"]
```



### ploidy

**type:** : `Int64`
>>>>>>> master

The ploidy of the samples

```julia
<<<<<<< HEAD
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
=======
2
```



### genotypes 

**type:** :  `::Dict` of `[loci] => Array{Tuple,1}`

The genotypes of the `loci`, as a `dictionary` of loci => genotypes. The loci are the dictionary keys `keys`, and the genotype `values` are an array of `tuples`, with each value corresponding to an allele. 

```
["locus_001"] => [(0,1),(0,0),(1,2)]
["locus_002"] => [(0,0),(1,1),(2,2)]
```

!!! important
    We use the **immutable** (cannot be changed) tuple type for genotypes of individuals because by the time you're using `PopGen.jl`, your data should already be filtered and screened. Hand-editing of genotype values is **strongly** discouraged, and we don't much like the idea of using this package to fudge your data that way.



### longitude

**type**: one-dimensional array of integers or floating point numbers `::Array{Union{Int64,Float64},1}`
>>>>>>> master

longitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```



<<<<<<< HEAD
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
=======
### latitude

**type**: one-dimensional array of integers or floating point numbers `::Array{Union{Int64,Float64},1}`

latitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```



## Viewing a PopObj

Given the volume of information that can be present in a `PopObj`, we recommend `summary()` to summarize the data rather than regurgitate everything on the screen. 
>>>>>>> master

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

   #Inds | Pop
   --------------
     21  |  1
     30  |  2
     28  |  3
     65  |  4
     28  |  5
     20  |  6
     20  |  7

Available fields: ind, popid, loci, ploidy, genotypes, longitude, latitude
```
<<<<<<< HEAD
=======



For context, this is what `gulfsharks` looks like without a semicolon at the end 🤮 :

![popobj_raw](img/popobj_raw.png)

We chose not to define a custom `Base.show` method to facilitate Atom/Juno visibility and drop-downs. 
>>>>>>> master
