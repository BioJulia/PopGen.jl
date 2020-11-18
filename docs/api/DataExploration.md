---
id: dataexploration
title: DataExploration.jl
sidebar_label: DataExploration.jl
---

### `missing`
```julia
missing(data::PopData; by::String = "sample")
```
Get missing genotype information in a `PopData`. Specify a mode of operation to return a DataFrame corresponding with that missing information.

**Modes**
- `"sample"` - returns a count and list of missing loci per individual (default)
- `"pop"` - returns a count of missing genotypes per population
- `"locus"` - returns a count of missing genotypes per locus
- `"full"` - returns a count of missing genotypes per locus per population

**Example**
```julia
missing(@gulfsharks, by = "pop")
```

-----

### `pairwise_identical`
    pairwise_identical(data::PopData)
Return a table of the percent of identical genotypes at each locus between all pairs of all individuals.

**Example**
```julia
julia> cats = @nancycats;
julia> pairwise_identical(cats)
27966×4 DataFrame
│ Row   │ sample_1 │ sample_2 │ identical │ n     │
│       │ String   │ String   │ Float64   │ Int64 │
├───────┼──────────┼──────────┼───────────┼───────┤
│ 1     │ N215     │ N216     │ 0.5       │ 8     │
│ 2     │ N215     │ N217     │ 0.25      │ 8     │
│ 3     │ N215     │ N218     │ 0.38      │ 8     │
│ 4     │ N215     │ N219     │ 0.38      │ 8     │
│ 5     │ N215     │ N220     │ 0.25      │ 8     │
│ 6     │ N215     │ N221     │ 0.5       │ 8     │
⋮
│ 27960 │ N296     │ N290     │ 0.0       │ 7     │
│ 27961 │ N297     │ N281     │ 0.14      │ 7     │
│ 27962 │ N297     │ N289     │ 0.43      │ 7     │
│ 27963 │ N297     │ N290     │ 0.29      │ 7     │
│ 27964 │ N281     │ N289     │ 0.25      │ 8     │
│ 27965 │ N281     │ N290     │ 0.43      │ 7     │
│ 27966 │ N289     │ N290     │ 0.14      │ 7     │
```

    pairwise_identical(data::PopData, sample_names::Vector{String})
Return a table of the percent of identical genotypes at each locus between all pairs of specified individuals.

**Example**
```julia
julia> cats = @nancycats;
julia> pairwise_identical(cats, samples(cats)[1:4])
6×4 DataFrame
│ Row │ sample_1 │ sample_2 │ identical │ n     │
│     │ String   │ String   │ Float64   │ Int64 │
├─────┼──────────┼──────────┼───────────┼───────┤
│ 1   │ N215     │ N216     │ 0.5       │ 8     │
│ 2   │ N215     │ N217     │ 0.25      │ 8     │
│ 3   │ N215     │ N218     │ 0.38      │ 8     │
│ 4   │ N216     │ N217     │ 0.12      │ 8     │
│ 5   │ N216     │ N218     │ 0.25      │ 8     │
│ 6   │ N217     │ N218     │ 0.0       │ 9     │
```