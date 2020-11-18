---
id: breeding_crosses
title: Simulate Breeding Crosses
sidebar_label: Breeding Crosses
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
import useBaseUrl from "@docusaurus/useBaseUrl";

:::note Requires PopGenSims.jl
To perfom simulations, you will need add and import the package `PopGenSims.jl` (available [here](https://github.com/pdimens/PopGenSims.jl)).
:::

If you need to simulate offspring genotypes given mating between two individuals, the `cross()` functions are available to simulate crosses and backcrosses.

**Currently, `PopGenSims.jl` can create crosses for:**
- haploids (ploidy = 1)
- diploids (ploidy = 2)
- tetraploids (ploidy = 4)
- hexaploids (ploidy = 6)
- octaploids (ploidy = 8)

## Perform a cross
```julia
cross(::PopData, parent1::String, parent2::String; n::Int, generation::String)
```
The cross function performs a simple parental cross from individuals `parent1` and `parent2` in the same PopData object. The parents are strings of the names of the parents in the PopData. The keyword argument `n` is the number of offspring to produce, and `generation` is a keyword argument for the population name to the assign the offspring (default: `"F1"`).

#### Example
```julia
julia> cats = @nancycats;

julia> f1 = cross(cats, "N111", "N107", n = 100000)
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100000
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent
```
Here is a look at the resulting `PopData`

<Tabs
  block={true}
  defaultValue="meta"
  values={[
    { label: 'meta', value: 'meta', },
    { label: 'loci', value: 'loci', },
  ]
}>
<TabItem value="meta">

There are two things that should jump out at you:
1. The `name` of offspring are prepended with `generation` and the `population` is the `generation`.
2. There is a never-before-seen `parents` column. This column exists for better record keeping of who has what parents if you are performing multiple crosses.

```
julia> f1.meta
100000×6 DataFrame
│ Row    │ name                │ ploidy │ population │ latitude │ longitude │ parents          │
│        │ String              │ Int64  │ String     │ Float32? │ Float32?  │ Tuple…           │
├────────┼─────────────────────┼────────┼────────────┼──────────┼───────────┼──────────────────┤
│ 1      │ F1_offspring_1      │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 2      │ F1_offspring_2      │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 3      │ F1_offspring_3      │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 4      │ F1_offspring_4      │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 5      │ F1_offspring_5      │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
⋮
│ 99995  │ F1_offspring_99995  │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 99996  │ F1_offspring_99996  │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 99997  │ F1_offspring_99997  │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 99998  │ F1_offspring_99998  │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 99999  │ F1_offspring_99999  │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
│ 100000 │ F1_offspring_100000 │ 2      │ F1         │ missing  │ missing   │ ("N111", "N107") │
```

</TabItem>
<TabItem value="loci">

```
julia> f1.loci
900000×4 DataFrame
│ Row    │ name                │ population │ locus  │ genotype   │
│        │ String              │ String     │ String │ Tuple…?    │
├────────┼─────────────────────┼────────────┼────────┼────────────┤
│ 1      │ F1_offspring_1      │ F1         │ fca8   │ (135, 135) │
│ 2      │ F1_offspring_1      │ F1         │ fca23  │ (146, 132) │
│ 3      │ F1_offspring_1      │ F1         │ fca43  │ (139, 145) │
│ 4      │ F1_offspring_1      │ F1         │ fca45  │ (132, 122) │
│ 5      │ F1_offspring_1      │ F1         │ fca77  │ (158, 150) │
⋮
│ 899995 │ F1_offspring_100000 │ F1         │ fca45  │ (122, 128) │
│ 899996 │ F1_offspring_100000 │ F1         │ fca77  │ (158, 150) │
│ 899997 │ F1_offspring_100000 │ F1         │ fca78  │ (142, 150) │
│ 899998 │ F1_offspring_100000 │ F1         │ fca90  │ (201, 199) │
│ 899999 │ F1_offspring_100000 │ F1         │ fca96  │ (113, 103) │
│ 900000 │ F1_offspring_100000 │ F1         │ fca37  │ (214, 208) │
```

</TabItem>
</Tabs>


## Perform a cross/backcross
```julia
cross(PopData => "Parent1Name", PopData => "Parent2Name", n::Int, generation::String)
```
This syntax uses the `Pair` notation of `PopData => "Parent"` to specify inputs. This method can be used for performing a cross like above, with the flexibility of parents allowed from two different `PopData` objects, which makes backcrosses possible. The keyword argument `n` is the number of offspring to produce, and `generation` is a keyword argument for the population name to the assign the offspring (default: `"F1"`).

#### Example
```julia
julia> f2_backcross = cross(cats => "N111", f1 => "F1_offspring_99", n = 100000, generation = "F2_manycats")
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100000
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent
```

And here you can see that `generation` was again prepended to each offspring `name`, along with assigned to the `population` for each.

```
julia> f2_backcross.meta
100000×6 DataFrame
│ Row    │ name                         │ ploidy │ population  │ latitude │ longitude │ parents                     │
│        │ String                       │ Int64  │ String      │ Float32? │ Float32?  │ Tuple{String,String}        │
├────────┼──────────────────────────────┼────────┼─────────────┼──────────┼───────────┼─────────────────────────────┤
│ 1      │ F2_manycats_offspring_1      │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 2      │ F2_manycats_offspring_2      │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 3      │ F2_manycats_offspring_3      │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 4      │ F2_manycats_offspring_4      │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 5      │ F2_manycats_offspring_5      │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
⋮
│ 99995  │ F2_manycats_offspring_99995  │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 99996  │ F2_manycats_offspring_99996  │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 99997  │ F2_manycats_offspring_99997  │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 99998  │ F2_manycats_offspring_99998  │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 99999  │ F2_manycats_offspring_99999  │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
│ 100000 │ F2_manycats_offspring_100000 │ 2      │ F2_manycats │ missing  │ missing   │ ("N111", "F1_offspring_99") │
```

:::caution
When crossing parents from different `PopData`, the parents must have the same loci. You will see error messages if they don't.
:::


## Merge results
The `PopData` generated from breeding crosses can be combined used `append` or `append!`

```julia
append(::PopData, ::PopData)
append!(::PopData, ::PopData)
```
These methods use outer joins and the `PopData` you are combining must have the same loci.

#### Example

```julia
# non mutating
crossed_sims = append(f1, f2_backcross)

# mutating
append!(f1, f2_backcross)
```