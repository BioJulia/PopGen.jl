---
id: read
title: Read.jl
sidebar_label: Read.jl
---

### `read_in`
```julia
read_in(infile::String; kwargs...)
```
Wraps `delimited()`, `genepop()`, `bcf()`, and `vcf()` to read a file in as a `PopData` object. File type is inferred from the file extension (case insensitive):
- delimited: `.csv` | `.tsv` | `.txt`
- genepop: `.gen` | `.genepop`
- variant call format: `.vcf` | `.bcf`

This function uses the same keyword arguments (and defaults) as the file importing functions it wraps; please see their respective docstrings in the Julia help console. (e.g. `?genepop`) for specific usage details. Use the alias function `file_import` interchangeably if you prefer the explicit name instead.

**Example**
```julia
read_in("cavernous_assfish.gen", digits = 3)
file_import("bos_tauros.csv", silent = true)
read_in("juglans_nigra.vcf")
```