Accessing any of these fields is done with a dot `.` accessor and can use the `[]`slice accessor, as per standard Julia convention:

```
julia> a = gulfsharks() ;

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
 
 julia> unique(a.popid)
 7-element Array{String,1}:
 "1"
 "2"
 "3"
 "4"
 "5"
 "6"
 "7"
```



See the "Accessing a PopObj" tutorial for some practice.

!!! failure "manual editing"

    Given the relationship of the ordered list of individuals, their population ID's, and the order of their genotypes in `.genotypes`, **NEVER** do any manual adding, removing, or sorting on `.ind`, `.popid`, or the arrays for any locus in `.genotypes`.