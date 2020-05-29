# Data exclusion

This section covers situations where one may want to remove samples or loci from `PopData`. By design, removal functions _do not_ alter your original `PopData`, so you always have a backup handy (but don't forget to assign a name to the new `PopData`). Since the exclusion functions don't alter the original `PopData`, they do not end with a bang `!`. 

::: details alias functions
For the exclusion commands on this page, we made it so the words `omit` and `remove` are interchangeable with `exclude`, meaning  `remove_samples` and `omit_loci` are the same exact functions as their `exclude` variants. This was done so you can use them interchangeably and you wont' need to remember the specific name to perform these basic actions. Maybe you just prefer the word `omit`. We're not here to judge.
:::

## Exclude samples

```julia
exclude_samples(data::PopObj, samp_id::Union{String, Vector{String}})
```

Returns a new `PopData` object without the sample or samples provided. Input can be a single sample, or an array of samples. Will output an entire `PopData`, so use a semicolon after the command if you don't want the entire object printed to your screen. Use `summary`  if you want to confirm that the samples were removed. This command will inform you if samples were not found in the data. 

:::: tabs card stretch
::: tab single individual
``` julia
julia> fewer_sharks = exclude_samples(sharks, "cc_001")
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 211
  Number of loci: 2213
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```
:::
::: tab multiple individuals
``` julia
julia> lots_fewer_sharks = remove_samples(sharks, ["cc_001", "cc_002", "cc_003"])
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 209
  Number of loci: 2213
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```
:::
::::
::: warning sample not found!
If removing a single sample and it is not found in the PopData, an error will be returned. However, if removing multiple samples, you will receive a notice above the PopData output indicating which individuals were not found, while still removing the ones that were present.
:::


## Exclude loci

```julia
exclude_loci(data::PopObj, loci::Union{String, Vector{String}})
```

Returns a new `PopData` object without the locus or loci provided. Input can be a single locus, or an array of loci, all as Strings. Will output an entire `PopData`, so use a semicolon after the command if you don't want the entire object printed to your screen. Use `summary`  if you want to confirm that the loci were removed. This command will inform you if loci were not found in the data.

:::: tabs card stretch
::: tab single locus
``` julia
julia> fewer_shark_loci = exclude_loci(sharks, "contig_475")
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 212
  Number of loci: 2212
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```
:::
::: tab multiple loci
``` julia
julia> lots_fewer_loci = remove_loci(sharks, ["contig_475", "contig_2784", "contig_8065"])
PopData Object
  Marker type: SNP
  Ploidy: 2
  Number of individuals: 212
  Number of loci: 2210
  Populations: 7
  Longitude: present with 0 missing
  Latitude: present with 0 missing
```
:::
::::

::: warning locus not found!
If removing a single locus and it is not found in the PopData, an error will be returned. However, if removing multiple loci, you will receive a notice above the PopData summary indicating which loci were not found, while removing the ones that were. If none of the loci specified were found, it will return an error.
:::