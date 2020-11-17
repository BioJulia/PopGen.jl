---
id: about
title: About
sidebar_label: About
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

PopGen.jl is an attempt to create a unified ecosystem for population genetics analyses leveraging the speed, power, fun(?), and community of the Julia language. We hope to implement common analyses (heterozygosity, kinship, FST, Tajima's D, etc.) in *sane*, user friendly ways, with consistent and comfortable syntax. 

:::note
The package is still in its infancy, so expect breaking changes to be more common than not between versions. We _think_ the core components are probably not going to change much, but cannot guarantee it. The package follows standard semantic versioning (i.e. breaking.feature.bugfix), so if the first number changes, we'll announce what's broken.
:::

## Goal

To be a comprehensive package for population genetics analyses and visualization that's fast and user friendly. This project is developed with a particular mantra: *Sanity, Sensibility, Accessibility*.

<Tabs
  block={true}
  defaultValue="sanity"
  values={[
    { label: 'Sanity', value: 'sanity', },
    { label: 'Sensibility', value: 'sensibility', },
    { label: 'Accessibility', value: 'accessibility', }
  ]
}>
<TabItem value="sanity">

Functions are written in a way such that their use is sane and natural. When possible (or sensible), we use full words for input variables or other components of input/output. The first argument for all functions is the input data without keywords. Always. If a function has a very obvious name, then it likely it won't need keyword arguments.

</TabItem>
<TabItem value="sensibility">

Functions need to be sensible, both in what they do and how they do it. This means they should include only the most relevant arguments and the most practical defaults. It also means the outputs need to be flexible enough to use with other Julia packages, such as `Query.jl`, or `Plots`.

</TabItem>
<TabItem value="accessibility">

Documentation is everything! Julia is a comparatively young language and we *want* you to be comfortable using PopGen.jl. That means investing **a lot** of time into writing thorough and visually pleasing documentation intended for users across a wide spectrum of proficiency. We also recognize that you (the reader) might not be very familiar with Julia, so we've written [a section](/getting_started/julia_primer.md) on clarifying some Julia concepts/conventions that will make reading this documentation a whole lot clearer.

</TabItem>
</Tabs>

## Why Julia (and not Python or R)?

### Speed and syntax

The speed can be comparable to C when coded using best practices. Also, has Python-like syntax with **optional** tabs. Not convinced? Check out our [comparison benchmarks](/getting_started/comparison.md). We are also fans of 1-indexing, but that's more of a perk.

### A modern language for modern problems

Julia has native support for: parallelization, distributed computing, GPU processing, and pipes! It also has robust machine learning packages (maybe for future work).

### Community & contribution

Julia's internals attempt to solve what's known as the "two language problem". That is, languages that are easy to write in (like Python, R, Ruby) are slow compared to languages that are more difficult to write in, which are fast (like C, C++, Fortran). For languages that are easier to write in, many of the commonly used packages and functions in those languages are written in another, faster language under the hood for performance reasons. On the whole, that's not really a problem for end-users, because things work and they are easy to use. **But**, it does become a problem when you want to investigate the code and implementation of a function. This means that, even as an R power-user, you are kind of helpless to investigate the implementation of something you are using in R if it's actually written in C++ under the hood. In a sense, it makes the publications of those methods less reproducible, because the users of it may be familiar with the language it's deployed in (like R), but not the language it's written in, (like C++). What if there are bugs?! What if the code implementation doesn't match the formulations in the publication?! **Yikes!**

So, if we write everything in Julia, and you use everything in Julia, users can themselves diagnose these things if they so choose. We want the community to have the power to make sure that all the math and logic checks out. This means anyone can contribute to the overall health and accuracy of this package. Yes, that means that we might be getting more Issues opened up (*ugh*), but a bug found is a lot better than a bug overlooked!

Like most Julia packages, PopGen.jl is written entirely in Julia, meaning the community using it need not learn another language if they wanted to contribute! Have you written a clever Julia function to calculate SAMOVA using `PopData`? Send us a pull request (please!)! Or [join the Slack group!](https://join.slack.com/t/popgenjl/shared_invite/zt-deam65n8-DuBs2z1oDtsbBuRplJW~Pg) 

#### Package manager

It's just such a delight to use. [Installing PopGen.jl](/getting_started/install.md) should be simple, consistent, and effortless thanks to the brilliant built-in package manager in Julia. 
