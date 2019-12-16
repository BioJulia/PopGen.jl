A little hands-on training will probably go a long way, so let's through some of the functions available in PopGen.jl with the included data. This tutorial will include both inputs and outputs so you can be confident what you're seeing in your Julia session is exactly what's supposed to happen. Sometimes the outputs can be a little lengthy, so they will be arranged in code "tabs", where the left-most tab is the input named after what it's accomplishing, and the right tab is the output of running the command. 



## Loading in the data

Let's keep things simple by loading in the nancycats data and calling it `ncats`

``` julia tab="load nancycats"
julia> ncats = nancycats() ; summary(ncats)
```

```  tab="output"
 Object of type PopObj:
 Marker type: Microsatellite
 Ploidy: 2
 Number of individuals: 237
 Number of loci: 9
 Longitude: none provided
 Latitude: none provided

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

Now that we have nancycats loaded in, we can use standard Julia accessor  conventions to view the elements within our PopObj.



## samples

### .samples

To view the entire `samples` dataframe

```julia tab=".samples"
julia> ncats.samples
```

```julia tab="output"
237×5 DataFrames.DataFrame
│ Row │ name   │ population    │ ploidy │ longitude │ latitude │
│     │ String │ Categorical…⍰ │ Int8   │ Any       │ Any      │
├─────┼────────┼───────────────┼────────┼───────────┼──────────┤
│ 1   │ N215   │ "1"           │ 2      │ missing   │ missing  │
│ 2   │ N216   │ "1"           │ 2      │ missing   │ missing  │
│ 3   │ N217   │ "1"           │ 2      │ missing   │ missing  │
│ 4   │ N218   │ "1"           │ 2      │ missing   │ missing  │
│ 5   │ N219   │ "1"           │ 2      │ missing   │ missing  │
│ 6   │ N220   │ "1"           │ 2      │ missing   │ missing  │
│ 7   │ N221   │ "1"           │ 2      │ missing   │ missing  │
│ 8   │ N222   │ "1"           │ 2      │ missing   │ missing  │
│ 9   │ N223   │ "1"           │ 2      │ missing   │ missing  │
│ 10  │ N224   │ "1"           │ 2      │ missing   │ missing  │
│ 11  │ N7     │ "2"           │ 2      │ missing   │ missing  │
│ 12  │ N141   │ "2"           │ 2      │ missing   │ missing  │
│ 13  │ N142   │ "2"           │ 2      │ missing   │ missing  │
│ 14  │ N143   │ "2"           │ 2      │ missing   │ missing  │
│ 15  │ N144   │ "2"           │ 2      │ missing   │ missing  │
│ 16  │ N145   │ "2"           │ 2      │ missing   │ missing  │
│ 17  │ N146   │ "2"           │ 2      │ missing   │ missing  │
⋮
│ 220 │ N258   │ "16"          │ 2      │ missing   │ missing  │
│ 221 │ N259   │ "16"          │ 2      │ missing   │ missing  │
│ 222 │ N260   │ "16"          │ 2      │ missing   │ missing  │
│ 223 │ N261   │ "16"          │ 2      │ missing   │ missing  │
│ 224 │ N262   │ "16"          │ 2      │ missing   │ missing  │
│ 225 │ N282   │ "17"          │ 2      │ missing   │ missing  │
│ 226 │ N283   │ "17"          │ 2      │ missing   │ missing  │
│ 227 │ N288   │ "17"          │ 2      │ missing   │ missing  │
│ 228 │ N291   │ "17"          │ 2      │ missing   │ missing  │
│ 229 │ N292   │ "17"          │ 2      │ missing   │ missing  │
│ 230 │ N293   │ "17"          │ 2      │ missing   │ missing  │
│ 231 │ N294   │ "17"          │ 2      │ missing   │ missing  │
│ 232 │ N295   │ "17"          │ 2      │ missing   │ missing  │
│ 233 │ N296   │ "17"          │ 2      │ missing   │ missing  │
│ 234 │ N297   │ "17"          │ 2      │ missing   │ missing  │
│ 235 │ N281   │ "17"          │ 2      │ missing   │ missing  │
│ 236 │ N289   │ "17"          │ 2      │ missing   │ missing  │
│ 237 │ N290   │ "17"          │ 2      │ missing   │ missing  │
```

### .name

This will access the names of the individuals as they appeared in the data.

``` julia tab=".ind"
julia> ncats.samples.name
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

### .population
This will access the names of the populations associated with each individual, in the same order as the  individuals.

``` julia tab=".popid"
julia> ncats.population
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

###  .ploidy

This shows you the ploidy of the data per individual.

``` julia tab=".ploidy"
julia> ncats.samples.ploidy
```

``` tab="output"
237-element Array{Int8,1}:
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 ⋮
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
 2
```

### .latitude

This accesses the latitude information of the PopObj. If there is none, like in the nancycats data, it returns an empty array.

```julia tab=".latitude"
julia> ncats.samples.latitude
```

```tab="output"
237-element Array{Any,1}:
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 ⋮      
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
```

### .longitude

This accesses the longitude information of the PopObj. Like before, if there is none, like in the nancycats data, it returns an array of missing values.

```julia tab=".longitude"
julia> ncats.samples.longitude
```

```tab="output"
237-element Array{Any,1}:
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 ⋮      
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
 missing
```

!!! info "Seeing some location info"
    The nancycats data has some weird coordinate system for information, so those data were omitted. If you want a proof of concept for `.longitude` and `.latitude`, load in `gulfsharks` and try it out. We'll use `hcat` (horizontal concatination) to horizontally bind the individual names, their latitude, and longitude. Later, you'll see that the `locations` command does this and a bit more.
    

``` julia tab="load gulfsharks"
julia> sharks = gulfsharks() ;    # semicolon just supresses printing output

julia> hcat(sharks.samples.name, sharks.samples.latitude, sharks.samples.longitude)
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

## loci

### .loci

This will show you the entire `loci` dataframe

```julia tab=".loci"
julia> ncats.loci
```

``` tab="output"
237×9 DataFrames.DataFrame. Omitted printing of 2 columns
│ Row │ fca23      │ fca37      │ fca43      │ fca45      │ fca77      │ fca78      │ fca8       │
│     │ Any        │ Any        │ Any        │ Any        │ Any        │ Any        │ Any        │
├─────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤
│ 1   │ (136, 146) │ (208, 208) │ (139, 139) │ (116, 120) │ (156, 156) │ (142, 148) │ missing    │
│ 2   │ (146, 146) │ (208, 208) │ (139, 145) │ (120, 126) │ (156, 156) │ (142, 148) │ missing    │
│ 3   │ (136, 146) │ (210, 210) │ (141, 141) │ (116, 116) │ (152, 156) │ (142, 142) │ (135, 143) │
│ 4   │ (138, 138) │ (208, 208) │ (139, 141) │ (116, 126) │ (150, 150) │ (142, 148) │ (133, 135) │
│ 5   │ (140, 146) │ (208, 208) │ (141, 145) │ (126, 126) │ (152, 152) │ (142, 148) │ (133, 135) │
│ 6   │ (136, 146) │ (208, 208) │ (145, 149) │ (120, 126) │ (150, 156) │ (148, 148) │ (135, 143) │
│ 7   │ (136, 146) │ (208, 208) │ (139, 145) │ (116, 126) │ (152, 152) │ (142, 148) │ (135, 135) │
│ 8   │ (136, 146) │ (208, 212) │ (135, 149) │ (120, 126) │ (154, 158) │ (142, 148) │ (135, 143) │
│ 9   │ (136, 146) │ (208, 212) │ (139, 139) │ (116, 126) │ (150, 160) │ (142, 142) │ (137, 143) │
│ 10  │ (132, 132) │ (208, 208) │ (141, 145) │ (120, 126) │ (150, 156) │ (148, 148) │ (135, 135) │
│ 11  │ (130, 136) │ (182, 182) │ (137, 145) │ (128, 128) │ (152, 152) │ (142, 150) │ (137, 141) │
│ 12  │ (130, 136) │ (182, 208) │ (135, 145) │ (126, 128) │ (144, 150) │ (140, 140) │ (129, 133) │
│ 13  │ (130, 130) │ (208, 208) │ (135, 145) │ (128, 130) │ (152, 156) │ (142, 142) │ (129, 133) │
│ 14  │ (130, 136) │ (182, 206) │ (135, 135) │ (128, 130) │ (156, 156) │ (142, 142) │ (133, 133) │
│ 15  │ (136, 136) │ (208, 208) │ (137, 137) │ (126, 130) │ (152, 152) │ (140, 142) │ (131, 135) │
│ 16  │ (136, 146) │ (182, 192) │ (135, 135) │ (128, 130) │ (144, 144) │ (142, 142) │ (129, 135) │
│ 17  │ (130, 144) │ (182, 192) │ (133, 133) │ (126, 126) │ (144, 144) │ (140, 140) │ (129, 133) │
⋮
│ 220 │ (136, 140) │ (208, 208) │ (139, 145) │ (122, 126) │ (144, 152) │ (142, 150) │ (137, 139) │
│ 221 │ (136, 140) │ (206, 208) │ (145, 149) │ (126, 126) │ (152, 152) │ (142, 142) │ (137, 143) │
│ 222 │ (136, 136) │ (210, 210) │ (145, 149) │ (122, 126) │ (152, 152) │ (150, 150) │ (135, 139) │
│ 223 │ (136, 136) │ (208, 208) │ (145, 145) │ (122, 126) │ (152, 156) │ (142, 142) │ (139, 143) │
│ 224 │ (136, 140) │ (208, 208) │ (149, 149) │ (120, 126) │ (152, 152) │ (142, 150) │ (135, 139) │
│ 225 │ (136, 138) │ (208, 208) │ (135, 139) │ missing    │ (150, 156) │ (142, 150) │ (133, 135) │
│ 226 │ (136, 136) │ (182, 182) │ (135, 139) │ missing    │ (146, 156) │ (142, 142) │ (133, 135) │
│ 227 │ (136, 136) │ (182, 208) │ (135, 139) │ missing    │ (150, 156) │ (142, 150) │ (133, 141) │
│ 228 │ (130, 146) │ (208, 208) │ (141, 141) │ missing    │ (148, 156) │ (142, 150) │ (133, 141) │
│ 229 │ (138, 138) │ (208, 208) │ (141, 145) │ missing    │ (148, 156) │ (142, 142) │ (123, 133) │
│ 230 │ (138, 138) │ (208, 208) │ (139, 139) │ missing    │ (150, 156) │ (142, 142) │ (123, 133) │
│ 231 │ (136, 146) │ (208, 208) │ (139, 139) │ missing    │ (150, 150) │ (142, 148) │ (133, 141) │
│ 232 │ (130, 136) │ (208, 208) │ (139, 145) │ missing    │ (152, 158) │ (142, 142) │ (133, 141) │
│ 233 │ (136, 146) │ (208, 220) │ (139, 145) │ missing    │ (150, 158) │ (142, 148) │ (133, 141) │
│ 234 │ (130, 130) │ (208, 208) │ (135, 145) │ missing    │ (148, 156) │ (142, 142) │ (133, 143) │
│ 235 │ (136, 144) │ (208, 208) │ (143, 143) │ missing    │ (144, 150) │ (142, 150) │ (135, 141) │
│ 236 │ (130, 136) │ (208, 208) │ (135, 145) │ missing    │ (150, 150) │ (142, 142) │ (137, 143) │
│ 237 │ (130, 146) │ (208, 208) │ (135, 139) │ missing    │ (150, 156) │ (142, 150) │ (135, 141) │
```



### locus names

This will access the names of the loci as they appeared in the data.

```julia tab=".loci"
julia> names(ncats.loci)
```

```tab="output"
9-element Array{Symbol,1}:
 :fca23
 :fca37
 :fca43
 :fca45
 :fca77
 :fca78
 :fca8 
 :fca90
 :fca96
```

You'll likely immediately notice the colons, which might not be what you expected, and that the type is `Array{Symbol,1}`. This is because `names` pulls the column names from the `.loci` dataframe of a `PopObj`, which are actually `Symbol` and not `String`. You can just as easily convert them to a string by broadcasting `String` over `names` using a dot `.`. This conversion is only for the output and won't change anything in the `PopObj` (nor does it need changing!)

```julia tab="conversion"
julia> String.(names(ncats.loci))
```

```tab="output"
9-element Array{String,1}:
 "fca23"
 "fca37"
 "fca43"
 "fca45"
 "fca77"
 "fca78"
 "fca8" 
 "fca90"
 "fca96"
```

### view genotypes

This is the core of the PopObj type. Each colum is an array of tuples that have the genotypes of each individual in the order with which they appear in `samples`. The convenience here is that each column of the dataframe is named for the locus, therefore you access genotypes with `PopObj.loci.locusname`

``` julia tab="access fca8"
julia> ncats.loci.fca8
```

``` tab="output"
237-element Array{Any,1}:
 missing    
 missing    
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

``` julia tab="access fca23"
julia> ncats.loci.fca23
```

``` tab="output"
237-element Array{Any,1}:
 (136, 146)
 (146, 146)
 (136, 146)
 (138, 138)
 (140, 146)
 (136, 146)
 (136, 146)
 (136, 146)
 (136, 146)
 (132, 132)
 (130, 136)
 (130, 136)
 (130, 130)
 (130, 136)
 (136, 136)
 (136, 146)
 (130, 144)
 (138, 138)
 (136, 144)
 (130, 136)
 ⋮         
 (136, 136)
 (136, 140)
 (136, 140)
 (136, 136)
 (136, 136)
 (136, 140)
 (136, 138)
 (136, 136)
 (136, 136)
 (130, 146)
 (138, 138)
 (138, 138)
 (136, 146)
 (130, 136)
 (136, 146)
 (130, 130)
 (136, 144)
 (130, 136)
 (130, 146)
```

## Slices
### general
You can likewise use slices to access parts of these data. If you're migrating from R or Python, it's the simple bracket accessor you're already familiar with, used to pull out a range of values. Julia just calls them slices. 

Let's look at a slice of `.name`.

``` julia tab="slice .name"
julia> ncats.samples.name[1:6]

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
julia> names(ncats.loci)[3:end]
```

``` tab="output"
7-element Array{Symbol,1}:
 :fca43
 :fca45
 :fca77
 :fca78
 :fca8 
 :fca90
 :fca96
```

!!! info ":end"
     All things start at 1, so there is no need for a special word for it. On the other hand, objects can have unknown lengths or varied lengths as you work with them. In Julia, use the word `end` in a slice range to indicate you want it to go to the end,, regardless of length or known size. 

One more example withh genotypes:

``` julia tab="slice genotypes"
julia> ncats.loci.fca8[1:3]
```

``` tab="output"
3-element Array{Any,1}:
 missing    
 missing    
 (135, 143)
```


## Operating on accessors

These accessors follow the exact same format as the dot operator in Python, or the `$` operator in R, meaning these objects can be assigned to new variables, you can operate on them, iterate over them, etc. 

Here's a simple example to display the unique population ID's in your PopObj

``` julia tab="unique"
julia> unique(ncats.samples.population)
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
    Julia has native piping (much like BASH) which uses the syntax `|>` (pipe + greater-than). With a pipe, `unique(ncats.samples.population)` can also be rewritten as `ncats.samples.population |> unique` for the same result. It's often a matter of preference for which you consider more readable. Use the `Pipe` package for even more robust piping where you can specify which argument the pipe relates to!

## 🛑❌ What to avoid! ❌🛑
Given the relationships of the ordered list of individuals (`.name`) and the order of genotypes in `.loci`, **NEVER USE `sort`, `sort!`, or manually arrange/add/delete anything in either dataframes `loci` or `samples`!!!** There are included functions `remove_loci!` and `remove_inds!` that do that kind of thing. That being said, you can rename individuals and their popid's without issue. Just no manual moving or deleting!