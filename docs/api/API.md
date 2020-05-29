# API

This page contains the APIs, or **A**pplication **P**rogramming **I**nterface, which are the entirety of all the functions/commands available in PopGen.jl. Unlike other sections of these docs, this page is intended to be *technical* rather than a guide. Included here are the function definitions and their docstrings as they appear inside this package. Most of these functions are used under-the-hood and not exported, meaning that if you want to use them, you will need to invoke them with `PopGen.function`. For example, if you wanted to use `unique_alleles` (which is not exported), you can do so with `PopGen.unique_alleles()` . 



## AlleleFreq.jl

:::: tabs
::: tab AlleleFreq.jl/alleles
### `alleles`
```julia
alleles(locus::T) where T<:GenoArray
```
Return an array of all the non-missing alleles present in a locus.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/unique_alleles
### `unique_alleles`
```julia
unique_alleles(locus::T) where T<:GenoArray
```
Return an array of all the unique non-missing alleles of a locus.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/allele_freq
### `allele_freq`
```julia
allele_freq(locus::T) where T<:GenoArray
```
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/allele_freq
### `allele_freq`
```julia
allele_freq(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of allele frequencies
of that locus per population.
:::
::: tab example
```julia
cats = nancycats()
allele_freq(cats, "fca8")
allele_freq(cats, "fca8", population = true)
```
:::
::::

:::: tabs
::: tab AlleleFreq.jl/allele_feq_vec
### `allele_freq_vec`
```julia
allele_freq_vec(locus::T) where T<:GenoArray
```
Return a Vector of allele frequencies of a single locus in a `PopData` object. Similar to `allele_freq()`, except it returns only the frequencies, without the allele names, meaning they can be in any order. This can be useful for getting the expected genotype frequencies.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_count_observed
### `geno_count_observed`
```julia
geno_count_observed(locus::T) where T<:GenoArray
```
Return a `Dict` of genotype counts of a single locus in a
`PopData` object.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_count_expected
### `geno_count_expected`
```julia
geno_count_expected(locus::T) where T<:GenoArray
```
Return a `Dict` of the expected genotype counts of a single locus in a
`PopData` object. Expected counts are calculated as the product of observed allele frequencies multiplied by the number of non-missing genotypes.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_freq
### `geno_freq`
```julia
geno_freq(locus::T) where T<:GenotypeArray
`PopData` object.
```
Return a `Dict` of genotype frequencies of a single locus in a
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_freq
### `geno_freq`
```julia
geno_freq(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.
:::
::: tab example
```julia
cats = nancycats()
geno_freq(cats, "fca8")
geno_freq(cats, "fca8", population = true)
```
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_freq
### `geno_freq`
```julia
geno_freq_expected(locus::T) where T<:GenoArray
```
Return a `Dict` of the expected genotype frequencies of a single locus in a `PopData` object. Expected frequencies are calculated as the product of
observed allele frequencies.
:::
::::

:::: tabs
::: tab AlleleFreq.jl/geno_freq_expected
### `geno_freq_expected`
```julia
geno_freq_expected(data::PopData, locus::String; population::Bool = false)
```
Return a `Dict` of expected genotype frequencies of a single locus in a
`PopData` object. Use `population = true` to return a table of expected genotype frequencies of that locus per population.
:::
::: example
```
cats = nancycats()
geno_freq_expeced(cats, "fca8")
geno_freq_expected(cats, "fca8", population = true)
```
:::
::::

## Datasets.jl

:::: tabs
::: tab Datasets.jl/dataset
### `dataset`
```julia
dataset(::String)
```
Load an example dataset from either `"gulfsharks"` (SNP) or `"nancycats"` (microsatellite). Can also use `"sharks"` and `"cats"`
as shorthands. Use `?nancycats` and `?gulfsharks` to learn more about
these datasets.
:::
::: tab example
```
ncats = dataset("nancycats")
gsharks = dataset("sharks")
```
:::
::::

:::: tabs
::: tab Datasets.jl/nancycats
### `nancycats`
```julia
nancycats()
```
Returns a `PopObj` of corresponding "nancycats" dataset as featured in
the R package `adegenet`. This is microsatellite data of 9 loci in 237
individuals across 17 populations.
:::
::: tab example
```
ncats = nancycats()
```
:::
::::

:::: tabs
::: tab Datasets.jl/gulfsharks
### `gulfsharks`
```julia
gulfsharks()
```
Returns a `PopObj` corresponding the Blacknose shark dataset as used in
Dimens et al. 2019. This is a mid-sized SNP dataset of 2213 SNPs in 212
individuals, across 7 populations.
:::
::: tab example
```julia
sharks = gulfsharks()
```
:::
::::


## HardyWeinberg.jl
:::: tabs
::: tab HardyWeinberg.jl/locus_chi_sq
### `locus_chi_sq`
```julia
locus_chi_sq(locus::T) where T <: GenoArray
```
Calculate the chi square statistic and p-value for a locus
Returns a tuple with chi-square statistic, degrees of freedom, and p-value.
:::
::::

:::: tabs
::: tab HardyWeinberg.jl/hwe_test
### `hwe_test`
```julia
    hwe_test(data::PopData; by::String = "locus", correction = "none")
```
Calculate chi-squared test of HWE for each locus and returns observed and expected heterozygosity with chi-squared, degrees of freedom and p-values for each locus. Use `by = "population"` to perform this separately for each population (default: `by = "locus"`). Use `correction =` to specify a P-value adjustment method for multiple testing.

**correction methods (case insensitive)**
- `"bonferroni"` : Bonferroni adjustment
- `"holm"` : Holm adjustment
- `"hochberg"` : Hochberg adjustment
- `"bh"` : Benjamini-Hochberg adjustment
- `"by"` : Benjamini-Yekutieli adjustment
- `"bl"`  : Benjamini-Liu adjustment
- `"hommel"` : Hommel adjustment
- `"sidak"` : Šidák adjustment
- `"forwardstop"` or `"fs"` : Forward-Stop adjustment
- `"bc"` : Barber-Candès adjustment
:::
::: tab example
```julia
hwe_test(gulfsharks(), correction = "bh")
hwe_test(gulfsharks(), by = "population", correction = "bh")
```
:::
::::

## Heterozygosity.jl
:::: tabs
::: tab Heterozygosity.jl/ishom
### `ishom`
```julia
ishom(locus::T) where T <: GenoArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if it is, `false` if it isn't, and `missing` if it's `missing`. The vector version simply broadcasts the function over the elements.
:::
::::

:::: tabs
::: tab Heterozygosity.jl/ishet
### `ishet`
```julia
ishet(locus::T) where T <: GenoArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if it is, `false` if it isn't. The vector version simply broadcasts the function over the elements. Under the hood, this function is simply `!ishom`.
:::
::::

:::: tabs
::: tab Heterozygosity.jl/hetero_o
### `hetero_o`
```julia
hetero_o(data::T) where T <: GenoArray
```
Returns observed heterozygosity as a mean of the number of heterozygous genotypes, defined as genotypes returning `true` for `ishet()`. This is numerically feasible because `true` values are mathematically represented as `1`, whereas `false` are represented as `0`.
:::
::::

:::: tabs
::: tab Heterozygosity.jl/hetero_e
### `hetero_e`
```julia
hetero_e(allele_freqs::Vector{T}) where T <: GenoArray
```
Returns the expected heterozygosity of an array of genotypes, calculated as 1 - sum of the squared allele frequencies.
:::
::::

:::: tabs
::: tab Heterozygosity.jl/heterozygosity
### `heterozygosity`
```julia
heterozygosity(data::PopData, mode::String = "locus")
```
Calculate observed and expected heterozygosity in a `PopData` object. For loci, heterozygosity is calculated in the Nei fashion, such that heterozygosity is calculated as the average over heterozygosity per locus per population.

**Modes**
- `"locus"` or `"loci"` : heterozygosity per locus (default)
- `"sample"` or `"ind"` or `"individual"` : heterozygosity per individual/sample
- `"population"` or `"pop"` : heterozygosity per population
:::
::: tab example
```julia
heterozygosity(nancycats(), "population" )
```
:::
::::

:::: tabs
::: tab Heterozygosity.jl/het_sample
### `het_sample`
```julia
    het_sample(data::PopData, individual::String)
```
Calculate the observed heterozygosity for an individual in a `PopData` object. Returns an array of heterozygosity values.
:::
::::

## Manipulate.jl

:::: tabs
::: tab Manipulate.jl/locations
### `locations`
```julia
    locations(data::PopData)
```
View the longitude and latitude data in a `PopData` object. Returns a table derived from the PopData. Changes made to this table will not alter the source `PopData` object.

Use `locations!` to add spatial data to a `PopData` object.
:::
::::

:::: tabs
::: tab Manipulate.jl/locations!
### `locations!`
```julia
locations!(data::PopData; lat::Vector{Float64}, long::Vector{Float64})
locations!(data::PopData, lat::Vector{Union{Missing,T}}, long::Vector{Union{Missing,T}}) where T <: AbstractFloat
locations!(data::PopData, lat::Vector{T}, long::Vector{T}) where T <: AbstractFloat
```
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes **decimal degrees** as a `Vector` of any `AbstractFloat` with tolerance for `missing` values.

## Formatting requirements
- Decimal Degrees format: `-11.431`
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)

:::
::: tab example

```julia
ncats = nancycats() ;
x = rand(237) ; y = rand(237)
locations!(ncats, long = x, lat = y)
```
:::
::::

:::: tabs
::: tab Manipulate.jl/locations!

### `locations!`
```julia
locations!(data::PopData; long::Vector{String}, lat::Vector{String})
```
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes **decimal minutes** or **degrees minutes seconds** format as a `Vector` of `String`. Recommended to use `CSV.read` from `CSV.jl` to import your spatial coordinates from a text file.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)



**NOTE**

If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, latitude degrees, latitude minutes), then the easiest course of action would be to merge them into two vectors of strings
(one for longitude, one for latitude):

```
long_string = string.(lat_deg, " ", lat_min)
lat_string = string.(long_deg, " ", long_min)
```
and use these as inputs into `locations!`

:::
::: tab example

```julia
ncats = nancycats();
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
locations!(ncats, long = x, lat = y)
```
:::
::::

:::: tabs
::: tab Manipulate.jl/loci

### `loci`
```julia
loci(data::PopData)
```
Returns an array of strings of the loci names in a `PopData` object.
:::
::::

:::: tabs
::: tab Manipulate.jl/get_genotypes

### `get_genotypes`
```julia
get_genotypes(data::PopObj; samples::Union{String, Vector{String}}, loci::Union{String, Vector{String}})
```
Return the genotype(s) of one or more `samples` for one or more specific `loci` (both as keywords) in a `PopData` object.
:::
::: tab example
```julia
cats = nancycats();
get_genotype(cats, samples = "N115" , loci = "fca8")
get_genotypes(cats, samples = ["N1", "N2"] , loci = "fca8")
get_genotype(cats, samples = "N115" , loci = ["fca8", "fca37"])
get_genotype(cats, samples = ["N1", "N2"] , loci = ["fca8", "fca37"])
```
:::
::::

:::: tabs
::: tab Manipulate.jl/get_sample_genotypes
### `get_sample_genotypes`
```julia
get_sample_genotypes(data::PopData, sample::String)
```
Return all the genotypes of a specific sample in a `PopData` object. This is an extension for the internal function `get_genotypes`.
:::
::: tab example
```julia
cats = nancycats()
get_sample_genotypes(cats, "N115")
```
:::
::::

:::: tabs
::: tab Manipulate.jl/locus
### `locus`
```julia
locus(data::PopData, locus::Union{String, Symbol})
```
Convenience wrapper to return a vector of all the genotypes of a single locus 
:::
::: tab example
```julia
locus(gulfsharks(), "contig_475")
```
:::
::::

:::: tabs
::: tab Manipulate.jl/missing
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
:::
::: tab example
```julia
missing(gulfsharks(), by = "pop")
```
:::
::::

:::: tabs
::: tab Manipulate.jl/populations
```julia
populations(data::PopData; listall::Bool = false)
```
View unique population ID's and their counts in a `PopData`.

- `listall = true` displays all samples and their `population` instead (default = `false`)
:::
::::

:::: tabs
::: tab Manipulate.jl/populations!
### `populations!`
```julia
populations!(data::PopData, rename::Dict)
```
Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`
:::
::: tab example
```julia
potatopops = Dict("1" => "Idaho", "2" => "Russet")
populations!(potatoes, potatopops)
```
:::
::::

:::: tabs
::: tab Manipulate.jl/populations!
### `populations!`
```julia
populations!(data::PopData, rename::Vector{String})
```
Replace the current population names with a `Vector` of new unique population names in the order that they appear in the PopData.meta.
:::
::: tab  example
```julia
potatopops = ["Idaho", "Russet"]
populations!(potatoes, potatopops)
```
:::
::::

:::: tabs
::: tab Manipulate.jl/populations!
### `populations!`
```julia
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Completely reassign populations for each individual. Takes two vectors of strings as input: one of the sample names, and the other with their new corresponding population name.

:::
::: tab example

```julia
populations!(potatoes, ["potato_1", "potato_2"], ["north_russet", "south_russet"])
```
:::
::::

:::: tabs
::: tab Manipulate.jl/exclude_loci

### `exclude_loci`
```julia
exclude_loci(data::PopData, locus::String)
exclude_loci(data::PopData, loci::Vector{String})
```
Exclude selected loci from a `PopData` object. Returns a new `PopData` object,
leaving the original intact. Synonymous with `omit_loci` and `remove_loci`.
:::
::: tab examples
```julia
new_cats = exclude_loci(nancycats(), "fca8")
very_new_cats = exclude_loci(nancycats(), ["fca8", "fca23"])
```
:::
::::

:::: tabs
::: tab Manipulate.jl/exclude_samples
### `exclude_samples`
```julia
exclude_samples(data::PopData, samp_id::String)
exclude_samples(data::PopData, samp_id::Vector{String})
```
Exclude selected samples from a `PopData` object. Returns a new `PopData` object, leaving the original intact. Synonymous with `omit_samples` and `remove_samples`.
:::
::: tab examples
```julia
exclude_samples(nancycats, "N100")
exclude_samples(nancycats, ["N100", "N102", "N211"])
```
:::
::::

:::: tabs
::: tab Manipulate.jl/samples
```julia
samples(data::PopData)
```
View individual/sample names in a `PopData`
:::
::::


## SummaryInfo.jl
:::: tabs
::: tab SummaryInfo.jl/allele_table
### `alele_table`
```julia
allele_table(data::PopData)
```
Return a "tidy" IndexedTable of the loci, their alleles, and their alleles' frequencies.
:::
::::

:::: tabs
::: tab SummaryInfo.jl/allele_avg
### `allele_avg`
```julia
    allele_avg(data::PopData, rounding::Bool = true)
```
Returns a NamedTuple of the average number of alleles ('mean') and standard deviation (`stdev`) of a `PopData`. Use `rounding = false` to not round results. Default (`true`) roundsto 4 digits.
:::
::: tab example
```julia
allele_avg(nancycats(), rounding = false)
```
:::
::::

:::: tabs
::: tab SummaryInfo.jl/richness
### `richness`
```julia
richness(data::PopData)
```
Calculates various allelic richness and returns a table of per-locus allelic richness. Use `populations = true` to calculate richness by locus by population.
:::
::::

:::: tabs
::: tab SummaryInfo.jl/summary
### `summary`
```julia
summary(data::PopData)
```
Prints a summary of the information contained in a PopData object
:::
::::

## Types.jl
:::: tabs
::: tab Types.jl/PopObj
### `PopObj`
```Julia
PopObj
```
Generic AbstractType for use in PopGen.jl
:::
::::

:::: tabs
::: tab Types.jl/PopData
### `PopData`
```
PopData(meta::IndexedTable, loci::IndexedTable)
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into a PopObj,
and instead should use the provided genepop, csv, or vcf file importers.

**`meta` ::IndexedTable** individual/sample data with the columns:

- `name` ::String the individual/sample names
- `population` ::String population names
- `ploidy` ::Int8 ploidy in order of `ind`
- `longitude` ::Float64 longitude values
- `latitude` ::Float64 latitude values

**`loci` ::IndexedTable** Long-format table of sample genotype records

- `name` ::CategoricalString the individual/sample names
- `population`::CategoricalString population name
- `locus` ::CategoricalString of locus name
- `genotype` Tuple of Int8 or Int16 depending on SNP or microsatellite
:::
::::


:::: tabs
::: tab Types.jl/GenoType
### `GenoType`
```julia
Genotype::DataType
```
For convenience purposes, an alias for `NTuple{N, <:Signed} where N`, which is the type describing individual genotypes in PopData.
:::
::::

:::: tabs
::: tab Types.jl/GenoArray
### `GenoArray`
```julia
GenoArray::DataType
```
For convenience purposes, an alias for an `AbstractVector` of elements `Missing` and `Genotype`, which itself is of type `NTuple{N, <:Signed} where N`. The definition as an `AbstractVector` adds flexibility for `SubArray` cases.
:::
::::

## Delimited.jl
:::: tabs
::: tab Delimited.jl/delimited
### `delimited`
```julia
delimited(infile::String; delim::Union{Char,String,Regex} = "auto", digits::Int64 = 3, silent::Bool = false)
```
Load a delimited-type file into memory as a `PopData` object. *There should be no empty cells
in your file*

**Arguments**
- `infile` : path to file
- `delim` : delimiter characters. By default uses auto-parsing of `CSV.File`
- `digits` : number of digits denoting each allele (default: `3`)
- `diploid`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent`   : whether to print file information during import (default: `true`)

**File formatting:**
- First row is column names in this order:
    1. name
    2. population
    3. longitude
    4. latitude
    5. locus_1_name
    6. locus_2_name
    7. etc...

**Missing data - Genotypes**

Missing genotypes can be formatted as all-zeros `000000` or negative-nine `-9`

**Missing data - Locations**

If location data is missing for a sample (which is ok!), make sure the value is written
as *0*, otherwise there will be transcription errors!

**Formatting example:**
```
name,population,long,lat,Locus1,Locus2,Locus3
sierra_01,mountain,11.11,-22.22,001001,002002,001001
sierra_02,mountain,11.12,-22.21,001001,001001,001002
snbarb_03,coast,0,0,001001,001001,001002
snbarb_02,coast,11.14,-22.24,001001,001001,001001
snbarb_03,coast,11.15,0,001002,001001,001001
```

:::
::: tab example
```julia
lizardsCA = delimited("CA_lizards.csv", digits = 3);
```
:::
::::

## Genepop.jl
:::: tabs
::: tab Genepop.jl/genepop
### `genepop`
```julia
genepop(infile::String; kwargs...)
```
Load a Genepop format file into memory as a PopData object.
- `infile::String` : path to Genepop file

**Keyword Arguments**
- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `diploid::Bool`  : whether samples are diploid for parsing optimizations (default: `true`)
- `silent::Bool`   : whether to print file information during import (default: `true`)

**File must follow standard Genepop formatting**
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default `POP`) must delimit populations
- Sample name is immediately followed by a *comma*
- File is *tab or space delimted* (but not both!)

**Genepop file example:**
```
wasp_hive.gen: Wasp populations in New York \n
Locus1
Locus2
Locus3
pop
Oneida_01,  250230  564568  110100
Oneida_02,  252238  568558  100120
Oneida_03,  254230  564558  090100
pop
Newcomb_01, 254230  564558  080100
Newcomb_02, 000230  564558  090080
Newcomb_03, 254230  000000  090100
Newcomb_04, 254230  564000  090120
```
:::
::: tab example
```julia
waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "pop")
```
:::
::::

## ioUtils.jl
:::: tabs
::: tab ioUtils.jl/determine_marker
### `determine_marker`
```julia
determine_marker(infile::String, geno_parse::CSV.File{}, digits::Int)
```
Return either `Int8` or `Int16` depending on largest allelic value in all genotypes in the first 10 samples of an input file (or all the samples if less than 10 samples). If the largest allele is 11 or greater, the marker will be considered a Microsatellite and coded in `PopData` as `Int16`, and the opposite is true for SNPs. There's no specific reason 10 was chosen other than it being a reasonable buffer for edge cases since SNP data <= 4, and haplotyped data could be a bit higher. Even if the microsatellite markers are coded incorrectly, there will be zero impact to performance, and considering how few microsatellite markers are used in typical studies, the effect on in-memory size should be negligible (as compared to SNPs).
:::
::::

:::: tabs
::: tab ioUtils.jl/find_ploidy
### `find_ploidy`
```julia
find_ploidy(genotypes::T where T<:SubArray)
```
Used internally in the `genepop` and `delimited` file parsers to scan the genotypes of a sample and return the ploidy of the first non-missing locus.
:::
::::

:::: tabs
::: tab ioUtils.jl/phase
### `phase`
```julia
phase(loc::String, type::DataType, digit::Int)
```
Takes a String of numbers returns a typed locus appropriate for PopGen.jl as used in the `genepop` and `csv` file parsers. Use `type` to specify output type (`Int8` or `Int16`), and `digit` to specify the number of digits/characters used per allele in a locus.
:::
::: tab example
```julia
ph_locus = phase("128114", Int16, 3)
map(i -> phase(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
:::
::::

:::: tabs
::: tab ioUtils.jl/phase_dip
### `phase_dip`
```julia
phase_dip(loc::String, type::DataType, digit::Int)
phase_dip(loc::T, type::DataType, digit::Int) where T<:Signed
phase_dip(loc::Missing, type::DataType, digit::Int)
```
Diploid-optimized variants of `phase()` that use integer division to split the alleles of a locus into a tuple of `type` Type. Use `type` to specify output type (`Int8` or `Int16`), and `digit` to specify the number of digits/characters used per allele in a locus. The `Missing` method simply returns `Missing`.
:::
::: tab examples
```julia
ph_locus = phase("128114", Int16, 3)
map(i -> phase_dip(i, Int16, 3), ["112131", "211112", "001003", "516500"])
# or #
[phase_dip(i, Int8, 2) for i in ["0101", "0103", "0202", "0103"]]
```
:::
::::


## Read.jl
:::: tabs
::: tab Read.jl/read_in
### `read_in`
```julia
read_in(infile::String; kwargs...)
```
Wraps `delimited())`, `genepop()`, `bcf()`, and `vcf()` to read a file in as a `PopData` object. File type is inferred from the file extension (case insensitive):
- delimited: `.csv` | `.tsv` | `.txt`
- genepop: `.gen` | `.genepop`
- variant call format: `.vcf` | `.bcf`
This function uses the same keyword arguments (and defaults) as the file importing functions it wraps; please see their respective docstrings in the Julia help console. (e.g. `?genepop`) for specific usage details. Use the alias function `file_import` interchangeably if you prefer the explicit name instead.
:::
::: tab example
```julia
read_in("cavernous_assfish.gen", digits = 3)
file_import("bos_tauros.csv", silent = true)
read_in("juglans_nigra.vcf")
```
:::
::::


## VariantCall.jl
:::: tabs
::: tab Variantcall.jl/bcf
### `bcf`
```julia
bcf(infile::String)
````
Load a BCF file into memory as a PopData object. Population and [optional] location information need to be provided separately.
- `infile` : path to BCF file
:::
::::

:::: tabs
::: tab VariantCall.jl/vcf
### `vcf`
```julia
vcf(infile::String)
```
Load a VCF file into memory as a PoDataj object. Population and [optional] location information need to be provided separately.
- `infile` : path to VCF file

:::
::::

## Utils.jl
:::: tabs
::: tab Utils.jl/convert_coord
```julia
convert_coord(coordinate::String)
```
Takes non-decimal-degree format as a `String` and returns it as a decimal degree `Float32`. Can be broadcasted over an array of coordinate strings to convert them.
### Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
:::
::: tab example
```julia
julia> convert_coord("-41 31.52")
-41.5253f0

julia> convert_coord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```
:::
::::

:::: tabs
::: tab Utils.jl/multitest_missing
### `multitest_missing`
```julia
multitest_missing(pvals::multitest_missing(pvals::Vector{T}, correction::String) where T <: Union{Missing, <:AbstractFloat}
```
Modification to `MultipleTesting.adjust` to include `missing` values in the returned array. Missing values are first removed from the array, the appropriate correction made, then missing values are re-added to the array at their original positions. See MultipleTesting.jl docs for full more detailed information.
:::
::: tab example
``` julia
multitest_missing([0.1, 0.01, 0.005, 0.3], "bh")
```
:::
::::

:::: tabs
::: tab Utils.jl/nonmissing
```julia
nonmissing(vec::T) where T<:AbstractArray
```
Convenience function to count the number of non-`missing` values
in a vector.
:::
::::

:::: tabs
::: tab Utils.jl/reciprocal
```julia
reciprocal(num::T) where T <: Signed
```
Returns the reciprocal (1/number) of a number. Will return `0` when the number is `0` instead of returning `Inf`.
:::
::::