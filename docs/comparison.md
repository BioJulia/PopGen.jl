# Comparison to `adegenet` / `pegas`

There's a reason we started investing so many hours and so many new grey hairs into writing PopGen.jl when there was an existing ecosystem in R to perform these same tasks. Like we explain in the home page of these docs, we want a platform that is 

1. fast(er)
2. written in a single language
3. easy to use

So, we'd like to prove that Julia and PopGen.jl actually achieves that by showing a few benchmarks comparing PopGen.jl to `adegenet` and `pegas`, which along with `ape` are arguably the most commonly used and robust population genetic packages available. It's worth mentioning that we ourselves use and have published with these packages, and are tremendously grateful for the work invested in those packages. We love you guys and respect the hell out of you! Here are links to [adegenet](https://github.com/thibautjombart/adegenet), [pegas](https://academic.oup.com/bioinformatics/article/26/3/419/215731/), and [ape](https://cran.r-project.org/package=ape).  



## Benchmarks

To make this a practical comparison, we're going to use the `gulfsharks` data because it is considerably larger (212 samples x 2213 loci) than `nancycats` (237 x 9) and a bit more of a "stress test".  All benchmarks in R are performed using the `microbenchmark` package, and  `BenchmarkTools` are used for Julia. 

```r tab="loading R packages"
library(adegenet)
library(pegas)
library(microbenchmark)
```

``` julia tab="loading julia packages"
using BenchmarkTools, PopGen
```

As a note, the reported benchmarks are being performed on a 64-bit Manjaro Linux system on a nothing-special Huawei Matebook D with 8gigs of RAM and a 4-core AMD Ryzen5 mobile processor. **Note:** none of the Julia benchmarks, unless explicitly stated, are using parallel or GPU processing. 



### Loading in data

Since `gulfsharks` is shamelessly provided in PopGen.jl, we simply invoke the `gulfsharks()` command in Julia. If you would like to try this yourself in R, find the `gulfsharks.gen` file in the package repository under `/data/data/gulfsharks.gen`. It will print out the input filename several times, which is omitted below for clarity.

```julia tab="Julia"
julia> @btime x = gulfsharks() ; # hide the output
  2.098 s (14175767 allocations: 707.58 MiB)
```

This R benchmark with take a few minutes. Consider going to the kitchen to make some tea while you wait.

```r tab="R"
> microbenchmark(read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L, quiet = TRUE))
Unit: seconds
                                                                                expr
 read.genepop(file = "/home/pdimens/gulfsharks.gen", ncode = 3L,      quiet = FALSE)
      min       lq     mean   median       uq      max neval
 5.670637 6.218719 6.745065 6.387936 7.019667 9.173005   100
```

Comparing averages, PopGen.jl clocks in at `2.098s` versus adegenet's `6.745s` , so ~3x faster.

Julia  :rocket:   |    R  :snail:


### `PopObj` vs `genind` size

It was pretty tricky to come up with a sensible/efficient/convenient data structure for PopGen.jl, and the original attempt was a Julian variant to a `genind`, which itself is something known as an `S4 class object`. While two dataframes design might not seem like it took a lot of effort, we ultimately decided that the column-major style and available tools, combined with careful genotype Typing was a decent "middle-ground" of ease-of-use vs performance. Plus, we are suckers for consistent syntax, which `genind`'s don't have compared to standard R syntax (looking at you too, Tidyverse/ggplot!). *Anyway*, it's important to understand how much space your data will take up in memory (your RAM) when you load it in, especially since data's only getting bigger! Keep in mind that `gulfsharks()` in PopGen.jl also provides lat/long data, which _should_ inflate the size of the object somewhat compared to the `genind`, which we won't add any location data to.

```julia tab="Julia"
julia> Base.summarysize(x)
4899793
```

vs

```r tab="R"
> object.size(gen)
5331536 bytes
```

![clutches pearls](img/clutches pearls cactus.png)

How is that possible?! Well, it's all in the Typing of the genotypes. Each genotype for each locus is encoded as a `Tuple` of either `Int8` (if SNPs) or `Int16` (if msats) to absolutely minimize their footprint without further going into byte-level encoding (so you can still see human-readable alleles).

The data takes up ~`4.9mb` in memory as a `PopObj` versus the ~`5.3mb` of a `genind`, about 18% smaller.

Julia  :house_with_garden: ​   |    R  :european_castle:


### Chi-squared test for HWE

This is a classic popgen test and a relatively simple one. 

```julia tab="Julia"
julia> @btime hwe_test(x, correction = "bh") ;
  479.550 ms (1954042 allocations: 74.12 MiB)
```

The R benchmark will take a while again, so if you're following along, this would be a good time to reconnect with an old friend.

```r tab="R"
> microbenchmark(hw.test(gen, B = 0))
Unit: seconds
                expr      min       lq     mean   median       uq      max neval
 hw.test(gen, B = 0) 5.100298 5.564807 6.265948 5.878842 6.917006 8.815179   100
```

Comparing averages, PopGen.jl clocks in at `480ms` versus adegenet's `6.3s`, so ~13x faster.

Julia  :rocket: |  R  :snail:

 