
## About

PopGen.jl is an attempt to shift population genetics analyses away from the patchwork of available pop-gen packages present in the R and Python languages, and combine it with the speed, power, fun(?), and community of the Julia language. We hope to implement common analyses (heterozygosity, kinship, FST, Tajima's D, etc.) in *sane*, user friendly ways, with syntax used within the package being consistent with the rest of the Julia ecosystem.

## Goal

To be a comprehensive package for population genetics analyses and visualization that's fast and user friendly. This project is developed with a particular mantra: *Sanity, Sensibility, Accessibility*.

:::: tabs card true
::: tab Sanity
Functions are written in a way such that their use is sane and natural. When possible (or sensible), we use full words for input variables or other components of input/output. The first argument for all functions is the input data without keywords. Always. If a function has a very obvious name, then it likely it won't need keyword arguments.
:::
::::

:::: tabs card true
::: tab Sensibility
Functions need to be sensible, both in what they do and how they do it. This means they should include only the most relevant arguments and the most practical defaults. It also means the outputs need to be flexible enough to use with other Julia packages, such as `Query.jl`, or `Plots`.
:::
::::

:::: tabs card true
::: tab Accessibility
Documentation is everything! Julia is a comparatively young language and we *want* you to be comfortable using PopGen.jl. That means investing **a lot** of time into writing thorough documentation intended for users across a wide spectrum of proficiency. We also recognize that you (the reader) might not be very familiar with Julia, so we've written [a section](/guide/) on clarifying some Julia concepts/conventions that will make reading this documentation a whole lot clearer.
:::
::::

## Why Julia (and not Python or R)?

#### Speed and syntax

The speed is comparable to C when coded correctly. Also, has Python-like syntax with **optional** tabs. We are also fans of 1-indexing, but that's more of a perk. Not convinced? Check out our [comparison benchmarks](getting_started/comparison.md).

#### A modern language for modern problems

Julia has native support for: parallelization, distributed computing, GPU processing, and pipes! It also has robust machine learning packages (maybe for future work).

#### Community & contribution

Julia's internals solve what's know as the "two language problem". That is, languages that are easy to write in (e.g. Python, R, Ruby) are slow compared to languages that are more difficult to write in, which are fast (e.g. C, C++, Fortran). For languages that are easier to write in, many of the commonly used packages and functions in those languages are written in another, faster language under the hood for performance reasons. On the whole, that's not really a problem for end-users, because things work and they are easy to use. **But**, it does become a problem when you want to investigate the code and implementation of a function. This means that, even as an R power-user, you are kind of helpless to investigate the implementation of something you are using in R if it's actually written in C++ under the hood. In a sense, it makes the publications of those methods less reproducible, because the users of it may be familiar with the language it's deployed in (like R), but not the language it's written in, (like C++). What if there are bugs?! What if the code implementation doesn't match the formulations in the publication?! **Yikes!**

So, if we write everything in Julia, and you use everything in Julia, users can themselves diagnose these things if they so choose. This means users can contribute to the overall health and accuracy of this package. Yes, that means that we might be getting more Issues opened up (*ugh*), but a bug found is a lot better than a bug overlooked!

Like most Julia packages, PopGen.jl is written entirely in Julia, meaning the community using it need not learn another language if they wanted to contribute! Have you written a clever Julia function to calculate SAMOVA using `PopData`? Send us a pull request (please!)! Or [join the Slack!](community.md)

#### Package manager

Which, let's be honest, is such a delight to use. [Installing PopGen.jl](getting_started/install.md) should be simple, consistent, and effortless thanks to the brilliant built-in package manager in Julia. 