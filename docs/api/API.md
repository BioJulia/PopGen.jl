---
id: api
title: API
sidebar_label: API
---

These pages contains the APIs, or **A**pplication **P**rogramming **I**nterface, which are the entirety of all the functions/commands available in PopGen.jl. Unlike other sections of these docs, these pages are intended to be *technical* rather than a guide. Included here are the function definitions and their docstrings as they appear inside this package. Most of these functions are used under-the-hood and not exported, meaning that if you want to use them, you will need to invoke them with `PackageName.function`. For example, if you wanted to use `unique_alleles` (which is not exported), you can do so with `PopGen.unique_alleles()`.

Here's a helpful map of navigating the API:
1. functions beginning with underscores are never exported
2. a red dot (🔴) indicates the function is not exported by any package
3. a pineapple (🍍) indicates the function is exported by PopGenCore.jl
4. a blue dot (🔵) indicates the function is exported by PopGen.jl

(this notation seemed easier than remembering 3 colors)