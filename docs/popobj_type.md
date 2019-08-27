For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called a `PopObj` (pronounced "pop ob" with a silent j because it rolls of the tongue better). If you haven't already guessed, it's a combination of the words Population and Object. The struct is defined as:

```julia
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genotypes::Dict
    longitude::Array{Union{Int64,Float64},1}
    latitude::Array{Union{Int64,Float64},1}
end
```

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
```



### loci

**type**: 1-dimension array of strings `::Array{String,1}`

The name of the loci, as an array of strings

```julia
["locus_001", "locus_2","super-awesome-locus-3"]
```



### ploidy

**type:** : `Int64`

The ploidy of the samples

```julia
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

longitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```



### latitude

**type**: one-dimensional array of integers or floating point numbers `::Array{Union{Int64,Float64},1}`

latitude data of samples (decimal degrees)

```
[-11.12, 15.32, 11.02, -4]
```



## Viewing a PopObj

Given the volume of information that can be present in a `PopObj`, we defined a custom `Base.show(::PopObj)` to summarize the data rather than regurgitate everything on the screen. 

```
julia> a
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



For context, this is what `/test/testdata.gen` looks like without a custom `show` function 🤮 :

![popobj_raw](img/popobj_raw.png)



