PopGen.jl provides two datasets as examples, each with their own easy-to-remember function for retrieving those data. 



## nancycats

We include the familiar nancycats data, as featured in `adegenet`, for easy download and import into PopGen.jl as a `PopObj`. Nancycats is a microsatellite dataset.

To use those data, simply invoke `nancycats()` without any arguments.

```
nancy = nancycats()
Object of type PopObj:
No location data provided

Number of individuals: 237
["N215", "N220", "N7"] … ["N294", "N296", "N281"]

Number of loci: 9
["fca8", "fca23", "fca43"] … ["fca90", "fca96", "fca37"]

Ploidy: 2
Number of populations: 4

   #Inds | Pop
   --------------
     64  |  1
     54  |  2
     67  |  3
     52  |  4

Available fields: ind, popid, loci, ploidy, genotypes, longitude, latitude
```

## gulfsharks

We also include the dataset using in Dimens *et al.* 2019 "[A **genomic** assessment of movement and gene flow around the South Florida vicariance zone in the migratory coastal **blacknose shark**, *Carcharhinus acronotus*](https://link.springer.com/article/10.1007/s00227-019-3533-1)" since it was already on hand. Like `nancycats`, we provide a convenient function to download and load these data into PopGen.jl as a `PopObj`. Gulfsharks is a SNP dataset.

To use those data, simply invoke `gulfsharks()` without any arguments. 

```jullia
sharkdata = gulfsharks()
Object of type PopObj:
No location data provided

Number of individuals: 212
["cca_001", "cca_002", "cca_003"] … ["seg_029", "seg_030", "seg_031"]

Number of loci: 2213
["contig_35208", "contig_23109", "contig_4493"] … ["contig_19384", "contig_22368", "contig_2784"]

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

