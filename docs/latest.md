---
id: latest
title: What's New
sidebar_label: What's New
---
## v0.7
### v0.7.0
#### ⚠️⚠️ Breaking Changes
- all PopData functionality moved to separate package PopGenCore.jl
- PopGen.jl reexports functions from PopGenCore.jl for familiar functionality
- `.meta` and `.loci` have been renamed `.metadata` and `.genodata`
- `.metadata` is no longer a DataFrame and instead a new `PopDataInfo` type
- lattitude and longitude columns no longer mandatory and omitted in cases where not used

#### ⚡⚡ Improvements
- PopData can be indexed like a DataFrame and it will return a brand new PopData!
- `PopDataInfo` is self-updating (in most cases)
- preliminary plink .bed file importing (not writing, yet)
- `show` for PopData is now smaller and cleaner
- INFO text for data importing now elides abs paths longer than the terminal width
- VCF/BCF support no longer lazy loaded
- VCF/BCF uses VariantCallFormat.jl now (instead of GeneticVariations.jl)
- VCF/BCF uses different GZ library for decompression
- `try...catch` blocks used in file io for faster file reading and fewer lines of code
-

#### 🐛🐛 Bug fixes
- super slow structure io on larger files


## v.0.6
### v.0.6.5
**Summary of changes from 0.6.2-5**

#### ⚡⚡ Improvements
- Bumped compat for `DataFrames.jl` to 1.0
- VCF/BCF importing now naturally sorts the loci names
  - includes new `NaturalSort.jl` dep
- file import `INFO` text consolidated somewhat
- `PopData` show method information consolidated somehwat

#### 🐛🐛 Bug fixes
- Hudson fst works as expected
- `isbiallelic` returns correct answer when used on `PopData` object
- [internal] conditional functions moved to `Conditionals.jl` file
- `keep` and `keep!` are exported


### v.0.6.1
#### ✨✨ New Features
- Hudson pairwise FST & Permutation
  - adds the Hudson et al. 1992 method
- `isbiallelic`
  - adds boolean test if a `PopData` object has only biallelic loci
  - adds boolean test if a `GenoArray` is biallelic
- `drop_multiallelic`
  - mutating and non-mutating methods to remove non-biallelic loci from a `PopData` object

#### ⚡⚡ Improvements
- `drop_monomorphic` now uses the same logic as `drop_multiallelic`, which should make it faster and leaner


#### 🐛🐛 Bug fixes
- `vcf` and `bcf` kwarg `rename_loci` now consistent in functions and docstrings
- `generate_meta` now uses a comprehension rather than deprecated `map(fn, groupeddataframe)` method


----

## v.0.5
### v.0.5.2
#### ⚡⚡ Improvements
- a rewrite of nei and weir-cockerham fst methods to be matrix-based (faster!)
#### ✨✨ New Features
- fully implements permutation testing for both pairwise fst methods
- adds method for `avg_allele_freq` to accommodate new `pairwise_nei`
- extends `pairwise_fst` to include iterations keyword to activate permutation testing

### v.0.5.1
#### ✨✨ New features
- `pairwise_fst` is now available for Weir & Cockerham (1984) and Nei (1987) methods
  - check out the [benchmarks](/docs/getting_started/comparison)!
- added `skipinf`, `skipnan`, and `skipinfnan` methods (unexported) to `Utils.jl`
- dropped `safemean` because the skip___ methods are a lot faster and slimmer


### v.0.5.0
This release fixes a critical bug in all the file importing functions that returned nothing when dropping monomorphic loci. Other changes include
- Dropping `JLD2.jl` suport due to its version-to-version instability. Two fewer dependencies!
  - As a result, `datasets()` now reads nancycats and gulfsharks directly from their source data files
  - To maintain all of the information, gulfsharks reads from a delimited file rather than a genepop file
- `populations` type signature and behavior has been changed:
  - the default returns an array of the unique population names
  - the keyword `listall::Bool` has been replaced with `counts::Bool`, which now returns a dataframe of the number of samples per population

-----

## v.0.4
### v0.4.5
This release builds off of `0.4.3` and does a better job with the VCF loading logic. Along with that, `vcf` and `bcf` exist in the namespace before loading in `GeneticVariation.jl`, meaning you can always view the docstrings. These stripped-down methods in the namespace will give helpful errors to remind you to load in `GeneticVariation.jl` and/or `GZip.jl`.

### v0.4.3
This release fixes and simplifies the under-the-hood `allele_freq`, `geno_freq`, and `geno_count_xxx` functions. The are faster now, and they infer types, making the output have expected type behavior. 

#### Changes
- You no longer need to import both `GeneticVariations.jl` and `GZip.jl` to have the `vcf` and `bcf` functions work. The reason is that if your file isn't gzipped, then why load in an unnecessary library? Therefore, if your file is gzipped, then you'll need to load in `GZip.jl` too, otherwise you just need `GeneticVariation.jl`. :cool:
- `avg_allele_freq` now has a different method, where the second positional argument is `power`, which will raise the calculated frequencies to the given value (default = `1`). This simplifies having to do things like square the values of the resulting `Dict`.

### v0.4.0
This release adds a slew of relatedness estimators, which can be bootstrapped and are performed in parallel. Paired with release of `PopGenSims.jl v0.0.2`.

#### ⚠️⚠️ Breaking changes
- CategoricalArrays replaced with PooledArrays
- VCF/BCF now lazy load and require `GeneticVariations.jl` and `GZip.jl` separately

#### ✨✨ New features
- relatedness estimators (see blog for tutorial)
- internal functions:
  - `loci_dataframe`
  - `loci_matrix`
  - `nonmissings`
  - `pairwise_pairs`
- `pairwiseidentical()` to compare percent identical loci
- `phase()` method
- Structure/fastStructure file IO

#### ⚡⚡ Improvements
- some internal function locations moved around (housekeeping)
- `nancycats()` and `gulfsharks()` are being phased out in favor of `@nancycats` and `@gulfsharks`. (You will see deprecation warning)
- documentation (Docusaurus) upgrades
  - edit button now correctly works on blog posts
- B/VCF reader rewritten (see docs)
