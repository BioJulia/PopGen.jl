PopGen.jl provides two datasets as examples, each with their own easy-to-remember function for retrieving those data. 



## nancycats

We include the familiar nancycats microsatellite data, as featured in `adegenet`, for easy importing into PopGen.jl as a `PopObj`.  

To use those data, simply invoke `nancycats()` without any arguments.

```
julia> ncats = nancycats() ; summary(ncats)
Object of type PopObj:
No location data provided

Number of individuals: 237
["N215", "N216", "N217"] … ["N281", "N289", "N290"]

Number of loci: 9
["fca8", "fca23", "fca43"] … ["fca90", "fca96", "fca37"]

Ploidy: 2
Number of populations: 17

Population names and counts:
17×2 DataFrames.DataFrame
│ Row │ population    │ count │
│     │ Categorical…⍰ │ Int32 │
├─────┼───────────────┼───────┤
│ 1   │ "1"           │ 10    │
│ 2   │ "2"           │ 22    │
│ 3   │ "3"           │ 12    │
│ 4   │ "4"           │ 23    │
│ 5   │ "5"           │ 15    │
│ 6   │ "6"           │ 11    │
│ 7   │ "7"           │ 14    │
│ 8   │ "8"           │ 10    │
│ 9   │ "9"           │ 9     │
│ 10  │ "10"          │ 11    │
│ 11  │ "11"          │ 20    │
│ 12  │ "12"          │ 14    │
│ 13  │ "13"          │ 13    │
│ 14  │ "14"          │ 17    │
│ 15  │ "15"          │ 11    │
│ 16  │ "16"          │ 12    │
│ 17  │ "17"          │ 13    │

Available .samples fields: .name, .population, .ploidy, .longitude, .latitude
```

The spatial coordinates provided for the dataset in `adegenet` are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted.  If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!

## gulfsharks

We also include the SNP dataset used in Dimens *et al.* 2019 "[A **genomic** assessment of movement and gene flow around the South Florida vicariance zone in the migratory coastal **blacknose shark**, *Carcharhinus acronotus*](https://link.springer.com/article/10.1007/s00227-019-3533-1)" since it was already on hand. Like `nancycats`, we provide a convenient function to load these data into PopGen.jl as a `PopObj`.

To use those data, simply invoke `gulfsharks()` without any arguments. 

```jullia
julia> sharks = gulfsharks() ; summary(sharks)
Object of type PopObj:

Longitude:
["-80.59928", "-80.59954", "-80.59958"] … ["-87.36617", "-85.71432", "-85.71432"]

Latitude:
["28.30624", "28.30787", "28.30234"] … ["30.05217", "29.82344", "29.82344"]


Number of individuals: 212
["cc_001", "cc_002", "cc_003"] … ["seg_029", "seg_030", "seg_031"]

Number of loci: 2213
["contig_35208", "contig_23109", "contig_4493"] … ["contig_19384", "contig_22368", "contig_2784"]

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

