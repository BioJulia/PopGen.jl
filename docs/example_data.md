PopGen.jl provides two datasets as examples, each with their own easy-to-remember function for retrieving those data. 



## nancycats

We include the familiar nancycats microsatellite data, as featured in `adegenet`, for easy importing into PopGen.jl as a `PopObj`.  

To use those data, simply invoke `nancycats()` without any arguments.

```
julia> ncats = nancycats() ; summary(ncats)
 Object of type PopObj
 Marker type: Microsatellite
 Ploidy: 2

 Number of individuals: 237
 Number of loci: 9
 Longitude: absent
 Latitude: absent

 Population names and counts:
17×2 DataFrame
│ Row │ population │ count │
│     │ Union…     │ Int64 │
├─────┼────────────┼───────┤
│ 1   │ P01        │ 10    │
│ 2   │ P02        │ 22    │
│ 3   │ P03        │ 12    │
│ 4   │ P04        │ 23    │
│ 5   │ P05        │ 15    │
│ 6   │ P06        │ 11    │
│ 7   │ P07        │ 14    │
│ 8   │ P08        │ 10    │
│ 9   │ P09        │ 9     │
│ 10  │ P10        │ 11    │
│ 11  │ P11        │ 20    │
│ 12  │ P12        │ 14    │
│ 13  │ P13        │ 13    │
│ 14  │ P14        │ 17    │
│ 15  │ P15        │ 11    │
│ 16  │ P16        │ 12    │
│ 17  │ P17        │ 13    │

```

The spatial coordinates provided for the dataset in `adegenet` are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted.  If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!

## gulfsharks

We also include the SNP dataset used in Dimens *et al.* 2019 "[A genomic assessment of movement and gene flow around the South Florida vicariance zone in the migratory coastal blacknose shark, *Carcharhinus acronotus*](https://link.springer.com/article/10.1007/s00227-019-3533-1)" since it was already on hand. Like `nancycats`, we provide a convenient function to load these data into PopGen.jl as a `PopObj`.

To use those data, simply invoke `gulfsharks()` without any arguments. 

```jullia
julia> sharks = gulfsharks() ; summary(sharks)
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

