# the PopObj type

## What is a PopObj

For the PopGen.jl package to be consistent, a standard flexible data structure needs to be defined. The solution is a custom type called a `PopObj` (pronounced "pop ob" with a silent j because it rolls of the tongue better).  The struct is defined as:

```julia
mutable struct PopObj
    ind::Array{String,1}
    popid::Array{Union{Int64,String},1}
    loci::Array{String,1}
    ploidy::Int64
    genotypes::Dict
#	hier			to be implemented later
    longitude::Array{Union{Int64,Float64},1}
    latitude::Array{Union{Int64,Float64},1}
end
```

- **ind** : the individual/sample names
- **popid** : the individual/sample population ID's (in the same order)
- **loci** : the name of the loci
- **ploidy** : the ploidy of the samples
- **genotypes** : the genotypes of the **loci**

- **longitude** : longitude data of samples (decimal degrees)
- **latitude** : latitude data of samples (decimal degrees)

Accessing any of these fields is done with a dot `.` accessor and can use the `[]`slice accessor, as per standard Julia convention:

```julia
julia> a = genepop("/test/testdata.gen", numpops = 7) ;

julia> a.ind[1:6]
6-element Array{String,1}:
 "cca_001"
 "cca_002"
 "cca_003"
 "cca_005"
 "cca_007"
 "cca_008"

julia> a.loci[1:6]
6-element Array{String,1}:
 "Contig_35208"
 "Contig_23109"
 "Contig_4493" 
 "Contig_10742"
 "Contig_14898"
 "Contig_8483" 

julia> a.ploidy
2
```



Given the volume of information that can be present in a `PopObj`, we defined a custom `Base.show(::PopObj)` to summarize the data rather than regurgitate everything on the screen. 

```julia
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



## Location Data

Notice the `No location data provided` text on the second line of output above. That text exists as a "heads up" rather than a warning because **location data is optional** for a `PopObj`. There are functions that use location information (e.g. `locations`and `plot_locations`), but most don't, so it's not a dealbreaker. If you add location information, displaying the `PopObj` again will show you output now including this information:

```julia
julia> a.latitude = rand(1:50, 212) ; a.longitude = rand(1:50, 212)	; # assign random location data

julia> a
Object of type PopObj:

Longitude:
["26", "7", "21"] … ["25", "26", "29"]

Latitude:
["10", "10", "40"] … ["3", "5", "46"]


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

