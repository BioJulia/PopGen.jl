---
id: add_info
title: Adding Information
sidebar_label: Adding Information
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

PopData's standard/default format includes information on:
- sample name
- sample population name
- sample ploidy
- sample geographical coordinates
- sample genotypes

But, sometimes you might want to add more information to the data structure. That's where the convenience function `add_meta!` comes in.

## `add_meta!`
This function has two methods, one for when the additional information you're adding is in the order with which your samples appear in `PopData.meta`, and another for when they don't. 

<Tabs
  block={true}
  defaultValue="s"
  values={[
    { label: 'same order', value: 's', },
    { label: 'different order', value: 'd', },
  ]
}>
<TabItem value="s">

```julia
add_meta!(popdata::PopData, metadata::T; name::String, loci::Bool = true, categorical::Bool = true) where T <: AbstractVector
```
Add an additional metadata information to a `PopData` object. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `PopData.meta`.

#### Arguments
- `popdata` : The `PopData` object to add information to
- `metadata` : A `Vector` with the metadata you wish to add to the `PopData`, in the same order as the names appear in `PopData.meta`

#### Keyword Arguments
- `name` : String of the name of this new column
- `loci` : Boolean of whether to also add this information to `PopData.loci` (default: `true`)
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `true`)

</TabItem>
<TabItem value = "d">

```julia
add_meta!(popdata::PopData, samples::Vector{String}, metadata::T; name::String, loci::Bool = true, categorical::Bool = true) where T <: AbstractVector
```
Add an additional metadata information to a `PopData` object. Mutates `PopData` in place. 
Takes a vector of sample names if the metadata is not in the same order as samples appear 
in `PopData.meta`.

#### Arguments
- `popdata` : The `PopData` object to add information to
- `sample` : A `Vector{String}` of sample names corresponding to the order of the provided `metadata` 
- `metadata` : A `Vector` with the metadata you wish to add to the `PopData`, in the same order as the names appear in `PopData.meta`

#### Keyword Arguments
- `name` : String of the name of this new column
- `loci` : Boolean of whether to also add this information to `PopData.loci` (default: `true`)
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `true`)

</TabItem>
</Tabs>

:::note futureproofing
As of yet, there are no features within PopGen.jl that require the use of `add_meta!`, but it is a great convenience function to have in our toolset for increasingly complicated things.
:::
