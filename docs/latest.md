---
id: latest
title: What's New
sidebar_label: What's New
---

## v0.4.0
This release adds a slew of relatedness estimators, which can be bootstrapped and are performed in parallel. Paired with release of `PopGenSims.jl v0.0.2`.

#### Breaking changes
- CategoricalArrays replaced with PooledArrays
- VCF/BCF now lazy load and require `GeneticVariations.jl` and `GZip.jl` separately

#### New features
- relatedness estimators (see blog for tutorial)
- internal functions:
  - `loci_dataframe`
  - `loci_matrix`
  - `nonmissings`
  - `pairwise_pairs`
- `pairwise_identical()` to compare percent identical loci
- `phase()` method
- Structure/fastStructure file IO


#### Changes
- some internal function locations moved around (housekeeping)
- `nancycats()` and `gulfsharks()` are being phased out in favor of `@nancycats` and `@gulfsharks`. (You will see deprecation warning)
- documentation (Docusaurus) upgrades
  - edit button now correctly works on blog posts
- B/VCF reader rewritten (see docs)