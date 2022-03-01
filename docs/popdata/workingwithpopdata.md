---
id: workingwithpopdata
title: Working with PopData
sidebar_label: Working with PopData
---

PopGen.jl includes basic commands to provide obvious methods to inspect and alter `PopData`. Using standard Julia conventions, only commands ending with a bang `!` are mutable, meaning they alter the input data. So, commands like `populations` will show you population information, whereas `populations!` will change that information in your `PopData`. The mutable commands here alter the data in your `PopData`, but not the source data (i.e. the files used to create the `PopData`). The "manipulation" commands were separated into smaller sections to make it less overwhelming, and using the `gulfsharks` data, you can explore each of the sections like a little tutorial. The sections don't follow any particular order, so feel free to jump around however you like. 

:::caution TLDR
End-users (vs developers) shouldn't access PopData fields directly
:::

In earlier versions of PopGen.jl, you were encouraged to directly access the internal fields of PopData. After careful consideration
and discussion with other users and developers, it's been decided that we should follow standard-ish convention and provide function
wrappers to view PopData fields and discourage direct access (unless you're a developer). This decision is intended to limit unintentional
errors, but also means a user has less to learn to get started.

To follow along like a tutorial, load the `@gulfsharks` or `@nancycats` data in if you haven't already:

```julia
julia> using PopGen

julia> sharks = @gulfsharks ;
```

## Topics
#### [Viewing Data](popdata/viewdata)
#### [Adding Information](popdata/addingdata)
#### [Data Exclusion](popdata/exclusion)
#### [Conditionals and Logic](popdata/conditionals)
#### [Population Data](popdata/populationdata)
#### [Location Data](popdata/locationdata)
#### [Data Exploration](popdata/dataexploration)
#### [Advanced Indexing](popdata/advancedindexing)
