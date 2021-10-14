---
id: exclusion
title: Data exclusion
sidebar_label: Data Exclusion
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

This section covers situations where one may want to selectively remove information from `PopData`. 

## Exclude certain things
```julia
exclude(data::PopData, kwargs...)
omit(data::PopData, kwargs...)
remove(data::PopData, kwargs...)
```
Returns a new `PopData` object without the criteria provided. Input can be a
single item or, or an array of items. You will be informed you if your criteria
were not found in the `PopData`. Using standard Julia conventions, `exclude()` will create a copy of your `PopData` excluding the specific criteria, whereas `exclude!()` will edit the input `PopData` in-place.

### Keyword Arguments
Everything gets converted to string, so `Symbol` works too if you want to cut down on keystrokes.
Integers work too if things are named `"1"`, `"2"`, etc.
- `locus::Union{String, Vector{String}}`: locus or loci you want to remove from the `PopData`
- `population::Union{String, Vector{String}}`: population(s) you want to remove from the `PopData`
- `name::Union{String, Vector{String}}`: sample(s) you want to remove from the `PopData`

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'samples', value: 's', },
    { label: 'loci', value: 'm', },
    { label: 'populations', value: 'p', },
    { label: 'combination', value: 'c', },
  ]
}>
<TabItem value="s">

``` julia
julia> fewer_sharks = exclude(sharks, name = "cc_001")
PopData Object
  Marker: SNP
  Ploidy: 2
  Samples: 211
  Loci: 2209
  Populations: 7
  Coordinates: present

julia> lots_fewer_sharks = remove(sharks, name = ["cc_001", "cc_002", "cc_003"])
PopData Object
  Marker: SNP
  Ploidy: 2
  Samples: 209
  Loci: 2209
  Populations: 7
  Coordinates: present
```

</TabItem>
<TabItem value="m">

``` julia
julia> fewer_shark_loci = exclude(sharks, locus = "contig_475")
PopData Object
  Marker: SNP
  Ploidy: 2
  Samples: 212
  Loci: 2208
  Populations: 7
  Coordinates: present

julia> lots_fewer_loci = remove(sharks, locus = ["contig_475", "contig_2784", "contig_8065"])
PopData Object
  Marker: SNP
  Ploidy: 2
  Samples: 212
  Loci: 2206
  Populations: 7
  Coordinates: present
```

</TabItem>
<TabItem value="p">

``` julia
julia> fewer_shark_pops = exclude(sharks, population = "Georgia")
PopData Object
  Markers: SNP
  Ploidy: 2
  Samples: 182
  Loci: 2209
  Populations: 6
  Coordinates: present

julia> lots_fewer_pops = remove(sharks, population = ["Florida Keys", "Mideast Gulf"])
PopData Object
  Markers: SNP
  Ploidy: 2
  Samples: 119
  Loci: 2209
  Populations: 5
  Coordinates: present
```

</TabItem>
<TabItem value="c">

``` julia
julia> tiny_shark = exclude(sharks, locus = "contig_475", name = ["cc_001", "neg_021", "mango_111"], population = ["Cape Canaveral", "kiwi"])
Notices:
  sample "mango_111" not found
  population "kiwi" not found

PopData Object
  Markers: SNP
  Ploidy: 2
  Samples: 190
  Loci: 2208
  Populations: 6
  Coordinates: present
```

</TabItem>
</Tabs>

The in-place variant `exclude!()` follows all the same syntax as `exclude()`, therefore all examples above would be identical for `exclude!()`.

:::note alias functions
The exclusion commands are interchangeable with `omit` and `remove`, both with and
without the bang (`!`). This was done so you can use the function comfortably without
needing to remember the specific name to perform it. Maybe you just prefer the word 
`omit` to `remove`. We're not here to judge. The examples below use any combination 
of `omit`, `remove`, and `exclude`.
:::

## Keep only certain things
```julia
keep(data::PopData, kwargs...)
```
Returns a new `PopData` object keeping only the occurrences of the specified keyword.
Unlike `exclude()`, only one keyword can be used at a time. All values are 
converted to `String` for filtering, so `Symbol` and numbers will also work.
### Keyword Arguments
- `locus::Union{String, Vector{String}}`: locus or loci you want to keep in the `PopData`
- `population::Union{String, Vector{String}}`: population(s) you want to keep in the `PopData`
- `name::Union{String, Vector{String}}`: sample(s) you want to keep in the `PopData`.

**Examples**
```
cats = @nancycats;
keep(cats, population = 1:5)
keep(cats, name = ["N100", "N102", "N211"])
keep(cats, locus = [:fca8, "fca37"])
```

## Remove monomorphic loci
While included in the file parsers by default, you may want to do this manually with
`drop_monomorphic`, which returns a new `PopData` object excluding any
monomorphic loci. You can use the mutable version `drop_monomorphic!` to
edit a `PopData` object in-place.
```julia
drop_monomorphic(::PopData)
drop_monomorphic!(::PopData)
```

## Remove multiallelic markers
If your data isn't biallelic, you may want to remove multi-allelic markers for
certain analyses (for example, the Hudson pairwise FST method requires
biallelic data). For that we have `drop_multiallelic`, which returns a new
`PopData` object, and the mutable version `drop_multiallelic!`, which edits a `PopData` object in-place.
```julia
drop_multiallelic(::PopData)
drop_multiallelic(::PopData)
```

**Example**
```
julia> sharks = @gulfsharks
PopData Object
  Markers: SNP
  Ploidy: 2
  Samples: 212
  Loci: 2209
  Populations: 7
  Coordinates: present


julia> drop_multiallelic(sharks)
[ Info: Removing 258 multialleic loci
PopData Object
  Markers: SNP
  Ploidy: 2
  Samples: 212
  Loci: 1951
  Populations: 7
  Coordinates: present
```