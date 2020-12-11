---
id: summaryinfo
title: SummaryInfo.jl
sidebar_label: SummaryInfo.jl
---

### `alele_table`
```julia
allele_table(data::PopData)
```
Return a "tidy" DataFrame of the loci, their alleles, and their alleles' frequencies.

----

### `allele_avg`
```julia
allele_avg(data::PopData, rounding::Bool = true)
```
Returns a NamedTuple of the average number of alleles ('mean') and standard deviation (`stdev`) of a `PopData`. Use `rounding = false` to not round results. Default (`true`) roundsto 4 digits.

**Example**
```julia
allele_avg(@nancycats, rounding = false)
```

----

### `richness`
```julia
richness(data::PopData)
```
Calculates various allelic richness and returns a table of per-locus allelic richness. Use `populations = true` to calculate richness by locus by population.

----

### `summary`
```julia
summary(data::PopData; by::String = "global")
```
Provides summary statistics for a `PopData` object. Use `by = "locus"` for
summary information by locus. Global values are given as unweighted means of 
the per-locus parameters.
#### Het_obs
observed heterozygosity given as:

1 - ∑ₖ ∑ᵢ Pₖᵢᵢ/np

where Pkii represents the proportion of homozygote `i` in sample `k` and `np` 
is the number of samples in that population
#### HT 
overall gene diversity given as:

1 - ∑ᵢ(p̄ᵢ² + (HS / (ñ × np)) - Het_obs / (2 × ñ × np))
where p̄ᵢ = ∑ₖpₖᵢ / np
#### HS        
within population gene diversity given as:

1 - ∑ᵢ(pᵢ² + HS / (ñ × np) - Het_obs / (2 × ñ × np))
where ñ = np / ∑ₖ(1/nₖ)
where p̄ᵢ² = ∑ₖ(pᵢₖ² / np)
#### DST       
amount of gene diversity among samples given as:

HT - HS
#### DST′      
amount of gene diversity among samples adjusted for sample size given as:

(np / (np-1)) × Dst
#### HT′       
overall gene diversity adjusted for sample size given as:

HS + DST′
#### FST       
proportion of the total genetic variance in subpopulations relative to the total genetic variance  given as:

DST / HT
#### FST′      
proportion of the total genetic variance in subpopulations relative to the total genetic variance, adjusted for 
heterozygosity given as:

DST′ / HT′
#### FIS       
proportion of the genetic variance in a locus relative to the genetic variance within subpopulations given as:

1 - (Het_obs / HS)
#### DEST      
population differentiation (Jost 2008) given as:

(np/(np-1)) × (Ht'-Hs) / (1-Hs)