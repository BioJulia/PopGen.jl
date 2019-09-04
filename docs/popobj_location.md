## Location Data

Notice the `No location data provided` text on the second line of output when you show a `PopObj`, like demonstrated on the previous page. That text exists as a "heads up" rather than a warning because **location data is optional** for a `PopObj`. There are functions that use location information (e.g. `locations`and `plot_locations`), but most don't, so it's not a dealbreaker. If you add location information, displaying the `PopObj` again will show you output now including this information:

```
julia> a = gulfsharks() ; summary(a)
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



!!! info
    While showing a `PopObj` represents location data as strings, they are actually coded as integers or floating point numbers. This was a design decision to remove visual clutter, however we can easily revert it back to displaying `Int` or `Float64` with type information if that's what users prefer.