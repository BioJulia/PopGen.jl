---
sidebarDepth: 2
---

# Comparison

There's a reason we started investing so many hours and so many new grey hairs into writing PopGen.jl when there was an existing ecosystem in R to perform these same tasks. Like we explain in the home page of these docs, we want a platform that is:

1. fast(er)
2. written in a single language
3. easy to use

So, we'd like to prove that Julia and PopGen.jl actually achieves that by showing a few benchmarks comparing PopGen.jl to `adegenet` and `pegas`, which along with `ape` are arguably the most commonly used and robust population genetic packages available. It's worth mentioning that we ourselves use and have published work incorporating these packages, and are tremendously grateful for the work invested in those packages. We appreciate those folks and have tremendous respect and envy for the work they continue to do! Here are links to [adegenet](https://github.com/thibautjombart/adegenet), [pegas](https://academic.oup.com/bioinformatics/article/26/3/419/215731/), and [ape](https://cran.r-project.org/package=ape).  



## Benchmarks

To make this a practical comparison, we're going to use the `gulfsharks` data because it is considerably larger (212 samples x 2213 loci) than `nancycats` (237 x 9) and a bit more of a "stress test".  All benchmarks in R are performed using the `microbenchmark` package, and  `BenchmarkTools` are used for Julia. 

:::: tabs card stretch
::: tab load R packages
```r
library(adegenet)
library(pegas)
library(microbenchmark)
```
:::
::: tab load Julia packages
``` julia
using BenchmarkTools, PopGen
```
:::
::::

As a note, the reported benchmarks are being performed on a 64-bit Manjaro Linux system on a nothing-special Lenovo Thinkbook 14S  with 8gigs of RAM and a 8th gen Intel i5 mobile processor. **Note:** all of the Julia benchmarks, unless explicitly stated, are performed single-threaded (i.e. not parallel, distributed, or GPU). 



### Loading in data

Since `gulfsharks` is shamelessly provided in PopGen.jl, we simply invoke the `gulfsharks()` command in Julia. If you would like to try this yourself in R, find the `gulfsharks.gen` file in the package repository under `/data/data/gulfsharks.gen`. It will print out the input filename several times, which is omitted below for clarity. Since the file importer now uses CSV.jl to read in the file, there are two steps of the genepop parser that are multithreaded. However, the majority of the data parsing (formatting the raw data into a correct PopObj structure) occurs using a single thread. This R benchmark will take a few minutes. Consider making some tea while you wait.

:::: tabs card stretch
::: tab Julia
```julia
julia> @btime x = gulfsharks() ; # hide the output
1.049 s (2472280 allocations: 166.30 MiB)
```
:::
::: tab R
```r
> microbenchmark(read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L, quiet = TRUE))
Unit: seconds
 read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L, quiet = FALSE)
      min       lq     mean   median       uq      max neval
 5.670637 6.218719 6.745065 6.387936 7.019667 9.173005   100
```
:::
::::

Comparing averages, PopGen.jl clocks in at `1.049s` versus adegenet's `6.745s` , so ~6.4x faster.

Julia  :rocket:   |    R  :snail:


### `PopData` vs `genind` size

It was pretty tricky to come up with a sensible/efficient/convenient data structure for PopGen.jl, and the original attempt was a Julian variant to a `genind`, which itself is something known as an `S4 class object`. While the two-IndexedTables design might not seem like it took a lot of effort, we ultimately decided that the column-major style and available tools, combined with careful genotype Typing was a decent "middle-ground" of ease-of-use vs performance. Plus, we are suckers for consistent syntax, which `genind`'s don't have compared to standard R syntax (looking at you too, Tidyverse/ggplot!). 

*Anyway*, it's important to understand how much space your data will take up in memory (your RAM) when you load it in, especially since data's only getting bigger! Keep in mind that `gulfsharks()` in PopGen.jl also provides lat/long data, which _should_ inflate the size of the object somewhat compared to the `genind`, which we won't add any location data to.
:::: tabs card stretch
::: tab Julia
```julia
julia> Base.summarysize(x)
1612428
#bytes
```
:::
::: tab R
```r
> object.size(gen)
5331536 bytes
```
:::
::::
![clutches pearls](/images/clutches_pearls_cactus.png)

What sorcery is this?! Well, it's all in the Typing of the genotypes. Each genotype for each locus is encoded as a `Tuple` of either `Int8` (if SNPs) or `Int16` (if msats) to absolutely minimize their footprint without further going into byte-level encoding (so you can still see human-readable alleles). An `Int8` is a signed integer that occupies 8bits of memory, whereas an `Int16` occupies 16bits (as compared to a standard `Int64`). 

The original file is `3.2mb`, and our `PopObj`takes up ~`1.6mb` in memory (half as big as the source file!) versus the ~`5.3mb` of a `genind`, which is ~1.5x larger than the source file and ~3.3x larger than the `PopData`. That's quite a big difference!

Julia  :house_with_garden: ​   |    R  :european_castle:


### Chi-squared test for HWE

This is a classic population genetics test and a relatively simple one. The R benchmark will take a while again, so if you're following along, this would be a good time to reconnect with an old friend.
:::: tabs card stretch
::: tab Julia
```julia
julia> @btime hwe_test(x, correction = "bh") ;
  392.527 ms (1599668 allocations: 57.20 MiB)
```
:::
::: tab R
```r
> microbenchmark(hw.test(gen, B = 0))
Unit: seconds
                expr      min       lq     mean   median       uq      max neval
 hw.test(gen, B = 0) 5.100298 5.564807 6.265948 5.878842 6.917006 8.815179   100
```
:::
::::
Comparing averages, PopGen.jl clocks in at ~`400ms` versus adegenet's `6.3s`, so ~15x faster.

Julia  :rocket:  |   R  :snail:

 
