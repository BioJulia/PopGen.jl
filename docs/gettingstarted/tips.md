---
id: tips
title: PopGen.jl tips
sidebar_label: PopGen.jl tips
---

Here are some useful tips to getting comfortable with PopGen.jl.

### PopGen vs PopGenCore
:::tip TLDR
PopGen is for analyses, PopGenCore is for everything else
:::

PopGen.jl is now intended exclusively for population genetic analyses, and PopGenCore.jl is the "core" package
where just about every utility lives relating to working with PopData that isn't an analysis. The recent split 
of PopGen and PopGenCore means PopGenCore is now a great standalone package for data viewing
and manipulation. If you don't need higher order analyses (which is what PopGen.jl provides), then PopGenCore.jl
should be enough. 


### Internal functions
:::tip TLDR
functions you aren't expected to use start with underscores `_`
:::

It's pretty common to have a series of user-facing functions in a package (the ones that get exported)
along with a series of unexported ones that are useful for development or are helper functions for the
exported ones. If you see a function that starts with an underscore, like `_adjacency_matrix`, then you
(as a user) aren't expected to know or worry about it, let alone use it. I mean, you totally can, but
not all of them have docstrings.

### Function names
:::tip TLDR
all user functions in PopGen.jl (not including PopGenCore.jl) are lowercase and use no underscores
:::

In an effort to be consistent with the Julia language and just consistent overall, user-facing functions 
are named with three main principals:
1. as verbose of a name as is reasonable
2. no underscores between words
3. all lowercase, always, unless it's a DataType

These decisions are because of years-long frustration with the swirling chaos that is R function names.
Where it is possible and reasonable to do so, function names are verbose and descriptive. For example,
if you wanted to perform a pairwise FST, the function is called `pairwisefst` -- **not** `fst`, `FST`, 
`pairwise_fst`, `pairwise_FST`, or `pairwiseFST`. This is an ongoing effort, so if you spot something
that doesn't fit this mold, please submit a Pull Request!

### Argument names
:::tip TLDR
all user functions in PopGen.jl have descriptive keyword argument names
:::

This too is driven by frustration with the shortened and/or truncated keyword arguments of just about
every R function. As if learning another language wasn't enough of a task, memorizing nonintuitive
keywords is a completely unneccessary stumbling block. To the best of our abilities, we try to name
keywords using _entire_ words, and hopefully words that are intuitive if you were just trying to guess
before calling up the docstring. 
