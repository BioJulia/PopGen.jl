# Installation

## Installing PopGen.jl

The package is currently unregistered while it's under active early development. However, installation is still simple!

Invoke the package manager with `]` in the REPL and use

```julia
add "https://github.com/pdimens/PopGen.jl"
```



If using Jupyter Notebooks or nteract, install the package in julia with

```julia
using Pkg
Pkg.add("https://github.com/pdimens/PopGen.jl") 
```



## Using PopGen

Like all Julia packages, to activate `PopGen`, simply run:

```julia
julia> using PopGen
```

Feel free to play around with the test data in `/test/testdata.gen` or add it to your workspace with the `nancycats` and `gulfsharks` commands.

!!! Note "Performance notes"
    If you're migrating to Julia from Python or R (or Matlab, etc.), you'll think Julia is slow and laggy because loading packages and running stuff has a noticeable wait time (10-40sec). However, if this is your first time in Julia, then it's worth mentioning that this lag is "compilation overhead". What this means is, Julia tries to pre-compile as much code as possible (into optimized machine code) when running something or loading a package. This lag exists **only the first time** you run something. Every subsequent run of a function, even with different parameters, will be **substantially** faster, and in most cases instant. If you want to test this yourself, try to run a line of code twice with `@time` before the function and compare the results. Here's an example:
    ``` 
    julia> @time using PopGen
    17.415902 seconds (19.88 M allocations: 1.022 GiB, 2.79% gc time)

    julia> @time using PopGen
    0.100233 seconds (64.07 k allocations: 3.123 MiB, 6.02% gc time)
    ```



