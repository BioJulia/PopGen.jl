---
id: install
title: Installation
sidebar_label: Installation
---

The package is currently unregistered while it's under active early development. However, installation is still simple!

:::: tabs card stretch
::: tab in the REPL or Juno
Invoke the package manager with `]` in the REPL and use

```julia
add "https://github.com/pdimens/PopGen.jl"
```

![install](/PopGen.jl/images/install_repl.gif)
:::
::: tab with Jupyter Notebooks or nteract

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/pdimens/PopGen.jl", rev="master"))
```

![install](/PopGen.jl/images/install_jupyter.gif)
:::
::::

## Using PopGen

Like all Julia packages, to activate PopGen.jl, simply run:

```julia
julia> using PopGen
```

Feel free to play around with the test data in `/data/data/` or add it to your workspace with the `nancycats` and `gulfsharks` commands.


::: warning Arch Linux users
If you compiled Julia from source, your PopGen.jl installation may fail due to incorrectly building `Arpack`, which is expected to be in one place, but the compilation puts in another. 

**Solutions**:

- install official Julia binaries from the AUR (`julia-bin`), which includes a correctly bundled `Arpack` (recommended)
- if Julia was compiled from source: install `julia-arpack` from the AUR and make sure to delete `~/.julia/packages/Arpack` if it exists. That *should* fix things, but sometimes it still acts up.
:::