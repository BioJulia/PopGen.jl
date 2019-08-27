# Editing PopObj elements

Following standard Julia convention, functions ending with `!`  are mutable, meaning they will edit the input. 

- To follow along, load in the provided test data in `/test/testdata.gen`

```julia
using PopGen
a = genepop("/test/testdata.gen", numpops = 7)
```

Keep in mind, your path to that file is probably different. **If using windows** make sure to change your backslashes "\" to forward slashes "/" 

## add location data

Location data can be added by directly accessing the fields `.longitude` and `.latitude` in your `PopObj`

```julia
a.longitude = rand(1:50, 212)   # creates 212 unique random numbers between 1 and 50
a.latitdue = rand(20:30, 212)	# creates 212 unique random numbers between 20 and 30
```



However, if your data is in decimal minutes rather than decimal degrees, use the `locations!` function to add it to the fields. This function will do a conversion for you.

###  adding decimal minutes data

```julia
locations!(x::PopObj; xloc::Array, yloc::Array)
```

Adds location data (longitude, latitude) to `PopObj`. Takes decimal degrees or decimal minutes format. **Must** use minus-sign instead of cardinal directions (i.e. 14 32.11W is **not** vaild). Location data must be in order of  individuals (`ind`). Replaces existing `PopObj` location data.

- Decimal Degrees : `-11.431`
- Decimal Minutes : `"-11 43.11"` (must use space and double-quotes)

If conversion is not necessary, can directly assign `PopObj.longitude` and `PopObj.latitude` as shown above.




## rename populations

```julia
popid!(x::PopObj; rename::Dict)
```
Rename the population ID's of `PopObj.popid`. Uses a `Dict` of `[popid] => replacement` to rename

Example:

```julia
# create a dictionary of name conversions
new_popnames = Dict(1 => "Cape Canaveral",
					2 => "Georgia",
					3 => "S Carolina",
    				4 => "FL Keys",
    				5 => "Mideast Gulf",
    				6 => "Northeast Gulf",
    				7 => "Southeast Gulf")

popid!(x, rename = new_popnames)
```



## remove loci

```julia
remove_loci!(x::PopObj, loci::Union{String, Array{String,1}})
```


Removes selected loci from a `PopObj`.

Examples:

```julia
remove_loci!(a, "contig_35208")

remove_loci!(a, ["contig_35208", "contig_23109", "contig_4493"])
```



## remove individuals

```julia
remove_inds!(x::PopObj, inds::Union{Array{String,1}})
```

Removes selected individuals from a `PopObj`.

Examples:

```julia
remove_inds!(a, "cca_001")

remove_inds(a, ["cca_001","cca_002", "cca_003"])
```