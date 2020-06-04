---
id: summaryinfo
title: SummaryInfo.jl
sidebar_label: SummaryInfo.jl
---

### `alele_table`
```julia
allele_table(data::PopData)
```
Return a "tidy" IndexedTable of the loci, their alleles, and their alleles' frequencies.

### `allele_avg`
```julia
    allele_avg(data::PopData, rounding::Bool = true)
```
Returns a NamedTuple of the average number of alleles ('mean') and standard deviation (`stdev`) of a `PopData`. Use `rounding = false` to not round results. Default (`true`) roundsto 4 digits.

**Example**
```julia
allele_avg(nancycats(), rounding = false)
```

### `richness`
```julia
richness(data::PopData)
```
Calculates various allelic richness and returns a table of per-locus allelic richness. Use `populations = true` to calculate richness by locus by population.
