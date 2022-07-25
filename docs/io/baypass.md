---
id: baypass
title: Baypass
sidebar_label: Baypass
---

[The Baypass software](http://www1.montpellier.inra.fr/CBGP/software/baypass/) is an increasingly common method to use to identify putative outlier loci in population datasets. The input file format is a matrix of loci (rows) x allele counts per population (columns). This format is not suitable for most other applications, so it cannot be read into PopGen.jl, but we offer a convenience function to write PopData into this format so you can use Baypass externally.

## baypass

```
baypass(data::PopData; filename::Union{String, Nothing} = nothing)
```
Convert a `PopData` object into a Baypass-format matrix. The required input format for the Baypass software
requires biallelic data. By default, it returns just the Baypass-format matrix; use the keyword argument `filename` to specify a file to write the matrix to.
This function **does not perform a Baypass analysis**, but instead creates the input matrix necessary for it.

The matrix specification is:
- rows = loci
    - each row is a different locus
- columns = allele counts per population
    - each pair of columns correspond to the alleles' counts (2 alleles, 2 columns) for a population
    - as a result, there should be 2 × n_populations columns
    - e.g. row 1, columns 1:2 are the allele counts for locus 1 in population 1

### Keyword arguments
- `filename`: a `String` of the name of the output file. If `nothing`, then this function just returns the Baypass input matrix without writing to a file. (default: `nothing`)

**Example**
```julia
julia> sharks = @gulfsharks;

julia> dropmultiallelic!(sharks)

julia> baypass(sharks, filename = "gulfsharks.baypass") ;
```