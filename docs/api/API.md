---
id: api
title: API
sidebar_label: API (start here)
---
import Icon from "@material-ui/core/Icon";

These pages contains the APIs, or **A**pplication **P**rogramming **I**nterface, which are the entirety of all the functions/commands available in PopGen.jl. Unlike other sections of these docs, these pages are intended to be *technical* rather than tutorial. Included here are the function definitions and their docstrings as they appear inside this package. Most of these functions are used under-the-hood and not exported, meaning that if you want to use them, you will need to invoke them with `PackageName.function`. For example, if you wanted to use `unique_alleles` (which is not exported), you can do so with `PopGen.unique_alleles()`. As a note, most of PopGenCore is reexported by PopGen.

Emoji will be used to indicate whether a function is exported or not and each API page
has a legend at the top for convenience. Here is some help making sense of it:
1. functions beginning with underscores are never exported
2. an eye with a cross (<Icon>visibility_off</Icon>) indicates the function is not exported by any package
3. a purple sqaure (🟪) indicates the function is exported by PopGenCore.jl
4. a blue dot (🔵) indicates the function is exported by PopGen.jl
5. a black dot (⚫) indicates the function is exported by PopGenSims.jl

