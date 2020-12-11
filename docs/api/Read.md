---
id: read
title: Read.jl
sidebar_label: Read.jl
---

### `read_from`
```julia
read_from(infile::String; kwargs...)
```
Wraps `delimited()`, `genepop()`, `bcf()`, and `vcf()` to read a file in as `PopData`. File type is
inferred from the file extension (case insensitive):

| File Format         | Extensions            | Docstring      |
| :------------------ | :-------------------- | :------------- |
| delimited           | `.csv`, `.txt`, `tsv` | `?delimited`   |
| genepop             | `.gen`, `.genepop`    | `?genepop`     |
| variant call format | `.vcf`, `.bcf`        | `?vcf`, `?bcf` |

This function uses the same keyword arguments (and defaults) as the file importing functions it wraps; please see their respective docstrings in the Julia help console. (e.g. `?genepop`) for specific usage details. Use the alias function `file_import` interchangeably if you prefer the explicit name instead.

**Example**
```julia
read_from("cavernous_assfish.gen", digits = 3)
file_import("bos_tauros.csv", silent = true)
read_from("juglans_nigra.vcf")
```

----

### `write_to`
```julia
write_to(data::PopData; filename::String, kwargs...)
```
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?popdata2filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the docstring to `popdata2genepop` with `?popdata2genepop`.

| File Format | Extensions             | Docstring            |
| :---------- | :--------------------- | :------------------- |
| genepop     | `.gen`, `.genepop`     | `?popdata2genepop`   |
| JLD2        | `.jld2`                | `?popdata2jld2`      |
| delimited   | `.csv`, `.txt`, `.tsv` | `?popdata2delimited` |

**Example**
```
cats = @nancycats;
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
write_to(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```

----

### `popdata2jld2`
```julia
popdata2jdl2(data::PopData; filename::String)
```
Write PopData to a `JLD2` file.
