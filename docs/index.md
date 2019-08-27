![logo](img/logo.png)



This is the documentation for the PopGen.jl package. If you're reading this, you're likely interested in doing some kind of population genetics analyses. Please read through the docs and try the tutorial out to get a feel for what PopGen.jl can do. 

## About

PopGen.jl is an attempt to shift population genetics analyses away from the patchwork of available pop-gen packages present in the R and Python languages, and combine it with the speed, power, fun(?), and community of the Julia language. We hope to implement common analyses (heterozygosity, kinship, FST, Tajima's D, etc.) in *sane*, user friendly ways, with syntax used within the package being consistent with the rest of the Julia ecosystem.



## Why Julia and not Python or R?

**Speed and syntax**

The speed is comparable to C when coded correctly. Also, has Python-like syntax with **optional** tabs. We are also fans of 1-indexing, but that's more of a perk.

**A modern language for modern problems**

Julia has native support for: parallelization, distributed computing, GPU processing, and pipes! It also has robust machine learning packages (maybe for future implementations).

**Community & contribution**

Julia's just-in-time type-safe and optimized compiling solves what's known as the "two language problem". That is, languages that are easy to write in are slow, and languages that are more difficult to write in are fast. Interpreted languages, like Python, R, and Perl are easier to write in, however many of the commonly used packages and functions in those languages are written in another, faster language under the hood for performance reasons. On the whole, that's not really a problem for end-users, because functions work and they are easy to use. **But**, it does become a problem when you want to investigate the code and implementation of a function. This means that, even as an R power-user, you are powerless to investigate the implementation of something you are using if it's actually written in C++. In a sense, it makes the publications of those methods less reproducible, because the users of it may be familiar with the language it's deployed in (like R), but not the language it's written in, (like C++). What if there are bugs?! What if the code implementation doesn't match the formulations in the publication?! **Yikes!**

And thus solving the two language problem, users can themselves diagnose code if they so choose. Yes, that means that we might be getting more Issues opened up, but a bug found is a lot better than a bug overlooked!

Like most Julia packages, PopGen.jl is written entirely in Julia, meaning the community using it need not learn another language if they wanted to contribute! Have you written a clever Julia function to calculate SAMOVA from a `PopObj`? Send us a pull request! (please!)

**Package manager**

Which is such a delight to use.

## Authors

Pavel Dimens, PhD Student @ U. Southern Mississippi

Jason Selwyn, PhD Candidate @ Texas A&M University - Corpus Christi

