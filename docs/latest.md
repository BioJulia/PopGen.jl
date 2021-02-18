---
id: latest
title: What's New
sidebar_label: What's New
---

## v0.4.5
This release builds off of `0.4.3` and does a better job with the VCF loading logic. Along with that, `vcf` and `bcf` exist in the namespace before loading in `GeneticVariation.jl`, meaning you can always view the docstrings. These stripped-down methods in the namespace will give helpful errors to remind you to load in `GeneticVariation.jl` and/or `GZip.jl`.

## v0.4.3
This release fixes and simplifies the under-the-hood `allele_freq`, `geno_freq`, and `geno_count_xxx` functions. The are faster now, and they infer types, making the output have expected type behavior. 

### Changes
- You no longer need to import both `GeneticVariations.jl` and `GZip.jl` to have the `vcf` and `bcf` functions work. The reason is that if your file isn't gzipped, then why load in an unnecessary library? Therefore, if your file is gzipped, then you'll need to load in `GZip.jl` too, otherwise you just need `GeneticVariation.jl`. :cool:
- `avg_allele_freq` now has a different method, where the second positional argument is `power`, which will raise the calculated frequencies to the given value (default = `1`). This simplifies having to do things like square the values of the resulting `Dict`.

-----

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
