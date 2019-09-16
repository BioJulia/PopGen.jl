A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs", where the left-most tab is the input named after what it's accomplishing, and the right tab is the output of running the command. 



## Loading in the data

Let's keep things simple by loading in the nancycats data and calling it `ncats`

``` julia tab="load nancycats"
julia> ncats = nancycats() ; summary(ncats)
```

```  tab="output"
Object of type PopObj:
No location data provided

Number of individuals: 237
["N215", "N216", "N217"] … ["N281", "N289", "N290"]

Number of loci: 9
["fca8", "fca23", "fca43"] … ["fca90", "fca96", "fca37"]

Ploidy: 2
Number of populations: 17

   #Inds | Pop
   --------------
     10  |  1
     22  |  2
     12  |  3
     23  |  4
     15  |  5
     11  |  6
     14  |  7
     10  |  8
     9   |  9
     11  |  10
     20  |  11
     14  |  12
     13  |  13
     17  |  14
     11  |  15
     12  |  16
     13  |  17

Available fields: ind, popid, loci, ploidy, genotypes, longitude, latitude
```



Note the last line of the output. These fields can be accessed by name using the dot `.` accessor.

## Trying out the standard accessors

Now that we have nancycats loaded in, we can use standard Julia accessor  conventions to view the elements within our PopObj.

### .ind
This will access the names of the individuals as they appeared in the data.

``` julia tab=".ind"
julia> ncats.ind
```

``` tab="output"
237-element Array{String,1}:
 "N215"
 "N216"
 "N217"
 "N218"
 "N219"
 "N220"
 "N221"
 "N222"
 "N223"
 "N224"
 "N7"  
 "N141"
 "N142"
 "N143"
 "N144"
 "N145"
 "N146"
 "N147"
 "N148"
 "N149"
 ⋮ 
 "N256"
 "N257"
 "N258"
 "N259"
 "N260"
 "N261"
 "N262"
 "N282"
 "N283"
 "N288"
 "N291"
 "N292"
 "N293"
 "N294"
 "N295"
 "N296"
 "N297"
 "N281"
 "N289"
 "N290"
```

### .popid
This will access the names of the populations associated with each individual, in the same order as the  individuals.

``` julia tab=".popid"
julia> ncats.ind
```

``` tab="output"
237-element Array{Union{Int64, String},1}:
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "1" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 "2" 
 ⋮   
 "16"
 "16"
 "16"
 "16"
 "16"
 "16"
 "16"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
 "17"
```

These ID's aren't super informative. Later, we'll change them using the `popid!` command.

### .loci
This will access the names of the loci as they appeared in the data.

``` julia tab=".loci"
julia> ncats.ind
```

``` tab="output"
9-element Array{String,1}:
 "fca8" 
 "fca23"
 "fca43"
 "fca45"
 "fca77"
 "fca78"
 "fca90"
 "fca96"
 "fca37"
```

###  .ploidy
This shows you the ploidy of the data, which is user-provided when loading in real data.

``` julia tab=".ploidy"
julia> ncats.ind
```

``` tab="output"
2
```

### .genotypes
This is the core of the PopObj type and it displays a `Dict` of loci corresponding to an array of tuples that have the genotypes of each individual in the order with which they appear in `.ind`. Due to the way `Dicts` parse (they do not have indices), viewing `.genotypes` will likely show you a different order of loci each time you display it. Worry not! Their values (the genotypes) are arrays, therefor indexed, and respect the order of `.ind`. Under the hood, PopGen.jl uses `.ind` and `.loci` to iterate over `.genotypes`, so ultimately, the order of loci is respected, even though there isn't any good analytical reason for it.

``` julia tab=".genotypes"
julia> ncats.genotypes
```

``` tab="output"
# the output scales to your REPL width

Dict{String,Array{Any,1}} with 9 entries:
  "fca77" => Any[(156, 156), (156, 156), (156, 152), (150, 150), (152, 152), (…
  "fca96" => Any[(113, 113), (113, 113), (113, 113), (91, 10, 5), (113, 113), …
  "fca23" => Any[(136, 146), (146, 146), (136, 146), (138, 138), (146, 140), (…
  "fca78" => Any[(142, 148), (142, 148), (142, 142), (142, 148), (142, 148), (…
  "fca8"  => Any[(0, 0), (0, 0), (135, 143), (135, 133), (135, 133), (135, 143…
  "fca45" => Any[(116, 120), (120, 126), (116, 116), (116, 126), (126, 126), (…
  "fca37" => Any[(208, 208), (208, 208), (210, 210), (208, 208), (208, 208), (…
  "fca43" => Any[(139, 139), (139, 145), (141, 141), (139, 141), (145, 141), (…
  "fca90" => Any[(199, 199), (199, 185), (197, 197), (199, 199), (199, 193), (…
```

Since `genotypes` is a `Dict`, you access it using standard Julia `Dict` conventions.  This includes the `keys` and `values` commands, or you can access the genotypes of specific loci with square brackets and the locus name in quotes:

``` julia tab="Dict at key"
julia> ncats.genotypes["fca8"]
```

``` tab="output"
237-element Array{Any,1}:
 (0, 0)    
 (0, 0)    
 (135, 143)
 (135, 133)
 (135, 133)
 (135, 143)
 (135, 135)
 (135, 143)
 (143, 137)
 (135, 135)
 (137, 141)
 (133, 129)
 (133, 129)
 (133, 133)
 (135, 131)
 (135, 129)
 (133, 129)
 (135, 129)
 (135, 135)
 (135, 131)
 ⋮         
 (139, 139)
 (137, 139)
 (137, 139)
 (143, 137)
 (135, 139)
 (143, 139)
 (135, 139)
 (135, 133)
 (135, 133)
 (133, 141)
 (133, 141)
 (133, 123)
 (133, 123)
 (133, 141)
 (133, 141)
 (133, 141)
 (143, 133)
 (135, 141)
 (143, 137)
 (135, 141)
```

###  .latitude
This accesses the latitude information of the PopObj. If there is none, like in the nancycats data, it returns an empty array.

``` julia tab=".latitude"
julia> ncats.latitude
```

``` tab="output"
0-element Array{Union{Float64, Int64},1}

```

###  .latitude
This accesses the longitude information of the PopObj. Like before, iIf there is none, like in the nancycats data, it returns an empty array.

``` julia tab=".longitude"
julia> ncats.longitude
```

``` tab="output"
0-element Array{Union{Float64, Int64},1}

```

!!! info "Seeing some location info"
    The nancycats data has some weird coordinate system for information, so those data were omitted. If you want a proof of concept for `.longitude` and `.latitude`, load in `gulfsharks` and try it out. We'll use `hcat` (horizontal concatination) to horizontally bind the individual names, their latitude, and longitude. Later, you'll see that the `locations` command does this and a bit more
    
    ``` julia tab="load gulfsharks"
    julia> sharks = gulfsharks() ;    # semicolon just supresses printing output
    
    julia> hcat(sharks.ind, sharks.latitude, sharks.longitude)
    ```
    
    ``` tab="output"
    212×3 Array{Any,2}:
    "cc_001"   28.3062  -80.5993
    "cc_002"   28.3079  -80.5995
    "cc_003"   28.3023  -80.5996
    "cc_005"   28.6123  -80.4225
    "cc_007"   27.8666  -80.3578
    "cc_008"   27.8666  -80.3579
    "cc_009"   27.8682  -80.3482
    "cc_010"   27.8711  -80.3482
    "cc_012"   28.4815  -80.4303
    "cc_013"   28.2421  -80.2494
    "cc_014"   28.2428  -80.2496
    "cc_015"   28.2378  -80.2496
    "cc_016"   28.2414  -80.2498
    "cc_017"   28.236   -80.2495
    "cc_018"   28.2364  -80.2494
    "cc_019"   28.3906  -80.4963
    "cc_020"   28.3902  -80.4968
    "cc_021"   28.3861  -80.4967
    "cc_022"   28.3869  -80.4966
    "cc_023"   28.3865  -80.496 
    ⋮                           
    "seg_001"  29.8901  -87.7189
    "seg_003"  30.1943  -88.0007
    "seg_009"  30.0021  -88.0493
    "seg_010"  30.0069  -88.049 
    "seg_011"  29.8362  -88.1675
    "seg_012"  29.5057  -88.0546
    "seg_014"  30.1428  -88.2974
    "seg_015"  30.2074  -88.36  
    "seg_016"  30.1151  -88.3922
    "seg_018"  29.8362  -88.168 
    "seg_021"  29.9466  -86.0399
    "seg_023"  29.9969  -85.6494
    "seg_024"  29.6966  -87.4403
    "seg_025"  30.064   -86.5374
    "seg_026"  30.0696  -86.5376
    "seg_027"  29.9065  -86.0905
    "seg_028"  30.0532  -87.3661
    "seg_029"  30.0522  -87.3662
    "seg_030"  29.8234  -85.7143
    "seg_031"  29.8234  -85.7143
    ```

## Slices
### general
You can likewise use slices to access parts of these data. If you're migrating from R or Python, it's the simple bracket accessor you're already familiar with, used to pull out a range of values. Julia just calls them slices. 

Let's look at a slice of `.ind`.

``` julia tab="slice .ind"
julia> ncats.ind[1:6]
```

``` tab="output"
6-element Array{String,1}:
 "N215"
 "N216"
 "N217"
 "N218"
 "N219"
 "N220"
```

Here's another example: 

``` julia tab="slice .loci"
julia> ncats.loci[3:end]
```

``` tab="output"
7-element Array{String,1}:
 "fca43"
 "fca45"
 "fca77"
 "fca78"
 "fca90"
 "fca96"
 "fca37"
```

!!! info ":end"
     All things start at 1, so there is no need for a special word for it. On the other hand, objects can have unknown lengths or varied lengths as you work with them. In Julia, use the word `end` in a slice range to indicate you want it to go to the end,, regardless of length or known size. 


### genotypes
We can do this for all the accessors, except `genotypes`, because it's a dictionary and has slightly different syntax. To use slices on `genotypes`, you will need to call the locus by name using the `Dict` syntax, then take a slice of that:

``` julia tab="slice genotypes"
julia> ncats.genotypes["fca8"][1:3]
```

``` tab="output"
3-element Array{Any,1}:
 (0, 0)    
 (0, 0)    
 (135, 143)
```





## Operating on accessors

These accessors follow the exact same format as the dot operator in Python, or the `$` operator in R, meaning these objects can be assigned to new variables, you can operate on them, iterate over them, etc. 

Here's a simple example to display the unique population ID's in your PopObj

``` julia tab="unique"
julia> unique(ncats.popid)
```

``` tab="output"
17-element Array{Union{Int64, String},1}:
 "1" 
 "2" 
 "3" 
 "4" 
 "5" 
 "6" 
 "7" 
 "8" 
 "9" 
 "10"
 "11"
 "12"
 "13"
 "14"
 "15"
 "16"
 "17"
```

!!! info "Pipes"
    Julia has native piping (much like BASH) which uses the syntax `|>` (pipe + greater-than). With a pipe, `unique(ncats.popid)` can also be rewritten as `ncats.popid |> unique` for the same result. It's often a matter of preference for which you consider more readable. Use the `Pipe` package for even more robust piping where you can specify which argument the pipe relates to!

## 🛑❌ What to avoid! ❌🛑
Given the relationships of the ordered list of individuals (`.ind`), their population ID's (`.popid`), and the order of genotypes in `.genotypes`, **NEVER USE `sort`, `sort!`, or manually arrange/add/delete anything in `.ind` or `popid`!!!** There are included functions `remove_loci!` and `remove_inds!` that do that kind of thing. That being said, you can rename individuals and their popid's without issue. Just no manual moving or deleting!

Let's look at what happens when you mananually delete an individual
```
julia> pop!(ncats.ind)
"N290"
```
Well, that seems benign, it removed the last individual from nancycats... Except it didn't, not fully. There is still a `popid` associated with that individual, and each of the 9 loci still have a genotype for it. If you had location data, it would still be there too. Now, the lengths of everything is messed up!
```
julia> length(ncats.ind)
236

julia> length(ncats.popid)
237

julia> length(ncats.genotypes["fca8"])
237
```

This is obviously a problem, and if you manually delete an individual from somewhere in the middle, then all the `popid`s and genotypes following it are shifted one, associated with the wrong individual. You would have to manually delete all the traces of this individual, which is annoying. Instead, did that work for you with the `remove_ind!` command!
