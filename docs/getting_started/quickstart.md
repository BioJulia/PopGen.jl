---
id: quickstart
title: Quick Start
sidebar_label: Quick Start
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## 1. Install PopGen

<Tabs
  block={true}
  defaultValue="repl"
  values={[
    { label: 'REPL/Juno', value: 'repl', },
    { label: 'Jupyter/nteract', value: 'jupyter', },
  ]
}>
<TabItem value="repl">

Invoke the package manager with `]` in the REPL and use

```julia
add PopGen
```

![install](/img/install_repl.gif)

</TabItem>
<TabItem value="jupyter">

Invoke `using Pkg`

```julia
using Pkg
Pkg.add("PopGen")
```

![install](/img/install_jupyter.gif)

</TabItem>
</Tabs>


## 2. Invoke the package
Simply run:

```julia
julia> using PopGen
```

## 3. Start playing around

To help get started, you can call `?PopGen` and be greeted with some information to help you get started:

```
julia> ?PopGen

        Quickstart for PopGen

Documentation: https://pdimens.github.io/PopGen.jl/

A few things things you can do to get started:

Load in data

- read_from(filename; kwargs...)
- genepop(infile; kwargs...)  or similar file-specific importer
- use available @gulfsharks or @nancycats datasets

Explore PopData

- populations(PopData) to view population names
- loci(PopData) to view locus names
- samples(PopData) to view sample names
- missing(PopData, by = ...) to view missing information

Manipulate PopData

- populations!(PopData, ...) to rename populations
- locations!(PopData, ...) to add geographical coordinates
- exclude!(PopData, kwargs...) to selectively remove data

Analyses

- richness(PopData) to calculate allelic richness
- relatedness(PopData) to calculate pairwise relatedness
- summary(PopData) to calculate F-statistics, heterozygosity, etc.
- hwe_test(PopData) to test for Hardy-Weinberg Equilibrium
```