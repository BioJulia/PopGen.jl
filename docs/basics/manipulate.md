---
id: manipulate
title: Start here
sidebar_label: Start here
---

PopGen.jl includes basic commands to provide obvious methods to inspect and alter `PopData`. Using standard Julia conventions, only commands ending with a bang `!` are mutable, meaning they alter the input data. So, commands like `populations` will show you population information, whereas `populations!` will change that information in your `PopData`. The mutable commands here alter the data in your `PopData`, but not the source data (i.e. the files used to create the `PopData`). Read over [Accessing parts of PopData](accessing) to become familiar with the components of `PopData`. 

The "manipulation" commands were separated into smaller sections to make it less overwhelming, and using the `gulfsharks` data, you can explore each of the sections like a little tutorial. The sections don't follow any particular order, so feel free to jump around however you like. 

To follow along like a tutorial, load the `gulfsharks` data in if you haven't already:

```julia
julia> using PopGen

julia> sharks = @gulfsharks ;
```

### [Accessing Elements](accessing)

### [Viewing Data](viewdata)

### [Sample and Locus Exclusion/Removal](exclusion)

### [Location and Population Information](populations)

