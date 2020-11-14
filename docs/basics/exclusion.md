---
id: exclusion
title: Data exclusion
sidebar_label: Data Exclusion
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

This section covers situations where one may want to selectively remove information from `PopData`. Using standard Julia conventions, `exclude()` will create a copy of your
`PopData` excluding the specific criteria, whereas `exclude!()` will edit the input
`PopData` in-place.

:::note alias functions
The exclusion commands are interchangeable with `omit` and `remove`, both with and
without the bang (`!`). This was done so you can use the function comfortably without
needing to remember the specific name to perform it. Maybe you just prefer the word `omit` to `remove`. We're not here to judge.

The examples below use any combination of `omit`, `remove`, and `exclude`, and are
similarly using varied keyword synonyms to demonstrate the flexibility of this
function's syntax, along with demonstrating that it doesn't matter whether the
keyword is singular or plural (like `locus` vs `loci`).
:::

## Exclude samples, loci, or entire populations

```julia
exclude(data::PopObj, kwargs...)
omit(data::PopObj, kwargs...)
remove(data::PopObj, kwargs...)
```
Returns a new `PopData` object without the sample or samples provided. Input can be a
single sample, or an array of samples. This command will inform you if your criteria
were not found in the `PopData`.

### Keyword Arguments
- `locus`: A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
    - The keyword `loci` also works.
- `population`: A `String` or `Vector{String}` of populations you want to remove from the `PopData`.
    - The keyword `populations` also works.
- `name`: A `String` or `Vector{String}` of samples you want to remove from the `PopData`.
    - The keywords `names`, `sample`, and `samples` also work.

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
julia> fewer_sharks = exclude(sharks, sample = "cc_001")
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
julia> fewer_shark_loci = exclude(sharks, loci = "contig_475")
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
julia> tiny_shark = exclude(sharks, loci = "contig_475", name = ["cc_001", "neg_021", "mango_111"], population = ["Cape Canaveral", "kiwi"])
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

:::note ______ not found!
If removing a single thing and it is not found in the PopData, an error will be returned. However, if attempting removing multiple things, you will receive a notice above the PopData output indicating which things were not found, while still removing the ones that were present.
:::

The in-place variant `exclude!()` follows all the same syntax as `exclude()`, therefore all examples above would be identical for `exclude!()`.

:::info Efficient Exclusion
1. Internally, there are two things that try to speed up data exclusion.
Exclusion occurs in the order of populations => samples => loci because this order removes the most volume of data per search so that each subsequent step has fewer data to process. Because of this order, you may get erroneous warnings if something isn't not found in the data if the previous step removed it.

2. If you are excluding more than half of the present data for any particular category, the function will instead just "keep the opposite". It's *much* faster to keep 53 loci than it is to remove 5400 of them (390ms vs 20s)!. 
:::