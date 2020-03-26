![logo](img/logo_slim.png)



This is the documentation for the PopGen.jl package. If you're here, you're likely interested in doing some kind of population genetics analyses. Please read through the docs and try the functions out with the test data to get a feel for what PopGen.jl can do. 

## About

PopGen.jl is an attempt to shift population genetics analyses away from the patchwork of available pop-gen packages present in the R and Python languages, and combine it with the speed, power, fun(?), and community of the Julia language. We hope to implement common analyses (heterozygosity, kinship, FST, Tajima's D, etc.) in *sane*, user friendly ways, with syntax used within the package being consistent with the rest of the Julia ecosystem.



## Goal

To be a comprehensive package for population genetics analyses and visualization that's fast and user friendly. This project is developed with a particular mantra: **Sanity, Sensibility, Accessibility**.

**Sanity**

Functions are written in a way such that their use is sane and natural. When possible (or sensible), we use full words for input variables or other components of input/output. The first argument for all functions is the input data without keywords. Always. If a function has a very obvious name, then it likely it won't need keyword arguments, such as `remove_loci!` , where it would be silly to have a keyword `loci=` since the name is explicit about what it does an the first argument is always the input data.

**Sensibility**

Functions need to be sensible, both in what they do and how they do it. This means they should include only the most relevant arguments, and the most sensible defaults. It also means the outputs need to be flexible enough to use with other Julia packages, such as `Query.jl`, or `Plots`.

**Accessibility**

Documentation is everything! We recognize Julia is a comparatively young language, and we know which languages and packages other population geneticists are using for their work. We *want* you to be comfortable using PopGen.jl, and that means investing a lot of time into writing thorough documentation intended for users (vs developers). We also recognize that you (the reader) might not be very familiar with Julia, or other non-R languages (which is ok!), so we've written a section on clarifying some Julia concepts/conventions that will make reading this documentation a whole lot clearer. It is by no means a replacement for sitting down and learning the Julia language a bit ([here is a great online book on that](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html)), but it should hopefully reduce the barrier of entry somewhat.

## Why Julia (and not Python or R)?

#### Speed and syntax

The speed is comparable to C when coded correctly. Also, has Python-like syntax with **optional** tabs. We are also fans of 1-indexing, but that's more of a perk. Not convinced? Check out our [comparison benchmarks](comparison.md).

#### A modern language for modern problems

Julia has native support for: parallelization, distributed computing, GPU processing, and pipes! It also has robust machine learning packages (maybe for future work).

#### Community & contribution

Julia's internals solve what's known as the "two language problem". That is, languages that are easy to write in (e.g. Python, R, Ruby) are slow compared to languages that are more difficult to write in, which are fast (e.g. C, C++, Java, Fortran). For languages that are easier to write in, many of the commonly used packages and functions in those languages are written in another, faster language under the hood for performance reasons. On the whole, that's not really a problem for end-users, because things work and they are easy to use. **But**, it does become a problem when you want to investigate the code and implementation of a function. This means that, even as an R power-user, you are kind of helpless to investigate the implementation of something you are using in R if it's actually written in C++ under the hood. In a sense, it makes the publications of those methods less reproducible, because the users of it may be familiar with the language it's deployed in (like R), but not the language it's written in (like C++). What if there are bugs?! What if the code implementation doesn't match the formulations in the publication?! **Yikes!**

So, if we write everything in Julia, and you use everything in Julia, users can themselves diagnose these things if they so choose. Yes, that means that we might be getting more Issues opened up, but a bug found is a lot better than a bug overlooked!

Like most Julia packages, PopGen.jl is written entirely in Julia, meaning the community using it need not learn another language if they wanted to contribute! Have you written a clever Julia function to calculate SAMOVA from a `PopObj`? Send us a pull request (please!)! Or [join the Slack!](community.md)

#### Package manager

Which, let's be honest, is such a delight to use. [Installing PopGen.jl](install.md) should be simple, consistent, and effortless thanks to the brilliant built-in package manager in Julia. 

## Authors

[![alt text](img/orcid.png)](https://orcid.org/0000-0003-3823-0373) [![alt text](img/twitter.png)](https://twitter.com/PVDimens) Pavel Dimens, PhD Student @ U. Southern Mississippi

[![alt text](img/orcid.png)](http://orcid.org/0000-0002-9100-217X) [![alt text](img/twitter.png)](https://twitter.com/JasonSelwyn) Jason Selwyn, PhD Candidate @ Texas A&M University - Corpus Christi 
