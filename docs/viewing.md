# Viewing elements in a PopObj

A handful of convenience functions exist for viewing data. 

- To follow along, load in the provided test data in `/test/testdata.gen`

```julia
using PopGen
a = genepop("/test/testdata.gen", numpops = 7)
```

Keep in mind, your path to that file is probably different. **If using windows** make sure to change your backslashes "\" to forward slashes "/" 



## individual/sample names

```julia
indnames(x::PopObj)
```

View individual/sample names in a `PopObj`. This is equivalent to `PopObj.ind`

Example:

```julia
indnames(a)
```



## loci

```julia
loci(x::PopObj, loci = nothing)
```

View the genotypes of all individuals for specific loci in a `PopObj`.
Default shows all genotypes for all individuals. Use `loci =` to specify a single
locus or array of loci to display.

Examples:

```julia
loci(a)

loci(a, "contig_35208")

loci(a, ["contig_35208", "contig_23109", "contig_4493"])
```



## locations

```julia
locations(x::PopObj)
```

View location data (`.longitude` and `.latitude`) in a `PopObj`

Example:

```julia
locations(a)
```

Use `locations!` to add spatial data to a `PopObj`



## population ID's

```julia
popid(x::PopObj; listall::Bool = false)
```

View unique population ID's in a `PopObj`. Default (`listall = false`) shows basic summary table of unique populations and the number of individuals in them.  

Use `listall = true` to display all individuals (`ind`) and their `popid` instead.

Examples:

```julia
popid(a)

popid(a, listall = true)
```

Using `popid!` to change the names of the populations



## genotypes

```julia
genotypes(x::PopObj; inds = nothing)
```

Show all the genotypes of specific individuals within a `PopObj`.  Default is to show all individuals so use `inds = ` to specify individuals.

- Names must be in quotes

Examples:

```julia
genotypes(a)

genotypes(a, inds = "cca_001")

genotypes(a, inds = ["cca_001", "cca_002", "cca_003"])
```





## find missing

```julia
missing(x::PopObj)
```



Identify and count missing loci in each individual of a `PopObj`. Returns a tuple of two 
`DataFrames`: loci per individual, number per loci.

Example:

```julia
missing(a)		# will display things

miss_ind, miss_loci = missing(a)
```

