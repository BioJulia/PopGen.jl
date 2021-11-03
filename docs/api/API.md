---
id: api
title: API
sidebar_label: API (start here)
---
import Icon from "@material-ui/core/Icon";

These pages contains the APIs, or **A**pplication **P**rogramming **I**nterface, which are the entirety of all the functions/commands available in PopGen.jl. Unlike other sections of these docs, these pages are intended to be *technical* rather than tutorial. Included here are the function definitions and their docstrings as they appear inside this package. Most of these functions are used under-the-hood and not exported, meaning that if you want to use them, you will need to invoke them with `PackageName.function`. For example, if you wanted to use `uniquealleles` (which is not exported), you can do so with `PopGenCore.uniquealleles()`. As a note, most things in PopGenCore are exported (for ease of development), and many parts of PopGenCore.jl are reexported by PopGen.jl.

Each API page features icons indicating whether a function is exported by that package. You will see two icons next to functions that PopGen.jl reexports from PopGenCore.jl.
Each page has a legend at the top for convenience. This is the icon system:
1. <Icon>visibility_off</Icon> indicates the function is not exported by any package
    - functions beginning with underscores always have <Icon>visibility_off</Icon>
3. <Icon>group_work</Icon> indicates the function is exported by PopGenCore.jl
4. <Icon>visibility</Icon> indicates the function is exported by the package the function appears in (PopGen.jl, PopGenSims.jl). 
    - this applies only to non-PopGenCore.jl packages

