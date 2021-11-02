---
id: locationdata
title: Location data
sidebar_label: Location data
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

Also known as geographic or coordinate data. The `sampleinfo` in the metadata can contain this information, which is not currently used in
any analyses present in PopGen.jl.

## View location data

```julia
locationdata(data::PopData)
```

View location (if present) data in a `PopData`,  returning a table the longitude and latitude information in the `metadata`. 

```julia
julia> locationdata(sharks)
212×2 SubDataFrame
 Row │ longitude  latitude 
     │ Float64    Float64  
─────┼─────────────────────
   1 │   28.3062  -80.5993
   2 │   28.3079  -80.5995
   3 │   28.3023  -80.5996
  ⋮  │     ⋮         ⋮
 210 │   30.0522  -87.3662
 211 │   29.8234  -85.7143
 212 │   29.8234  -85.7143
           206 rows omitted
```

## Add geographical coordinates
### decimal minutes
```julia
locationdata!(data::PopData; longitude::Vector{T}, latitude::Vector{T}) where T<:AbstractFloat
```
Location data can be added using one of the methods of `locations!`. As indicated by the bang `!`, your `PopData` will be edited in place, and there will be no return output. If your data is in anything other than Decimal-Degrees format, this function will convert your long/lat into Decimal Degrees. To import those data into Julia, you'll likely want to use the wonderful `CSV.jl` package first. The functions accept keywords `longitude` and `latitude`, or can be used without them so long as the vectors are input in that order. 

This method is pretty straightforward and tolerates vectors with `missing` data.
#### formatting requirements
- Coordinates must be decimal-minutes as either `Float32` or `Float64` (e.g. `-21.321`)

```julia
# generate some fake location data
julia> long = rand(212) .* 10 ; lat = rand(212) .* -10

julia> locationdata!(sharks, long, lat)
```


### other formats
```julia
locationdata!(data::PopData; longitude::Vector{String}, latitude::Vector{String})
```
#### formatting requirements

- Coordinates as `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)

If not already in decimal-minutes format, it would likely be most convenient if you imported your coordinate data as vectors of strings, which would look something like this:

```julia
longitude = ["-43 54 11", "22 23 11N"]
latitude = ["11 44 31", "-25 41 13"]
```

:::caution Missing values
This method tolerates `missing` values, but you will need to `replace!` instances of `missing` with the string `"missing"`.
:::