---
id: read
title: ReadWrite.jl
sidebar_label: ReadWrite.jl
---
## PopGenCore.jl/src/io/ReadWrite.jl
❗ => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### `read_from`
```julia
read_from(infile::String; kwargs...)
```
Wraps `delimited()`, `genepop()`, `structure()`, `bcf()`, and `vcf()` to read a file in as `PopData`. File type is
inferred from the file extension (case insensitive):

| File Format         | Extensions            | Docstring      |
| :------------------ | :-------------------- | :------------- |
| delimited           | `.csv`, `.txt`, `tsv` | `?delimited`   |
| genepop             | `.gen`, `.genepop`    | `?genepop`     |
| Structure/fastStructure | `.str`, `.structure` | `?structure` |
| variant call format | `.vcf`, `.bcf`        | `?vcf`, `?bcf` |

This function uses the same keyword arguments (and defaults) as the file importing functions it wraps; please see their respective docstrings in the Julia help console. (e.g. `?genepop`) for specific usage details. Use the alias function `file_import` interchangeably if you prefer the explicit name instead.

**Example**
```julia
read_from("cavernous_assfish.gen", digits = 3)
file_import("bos_tauros.csv", silent = true)
read_from("juglans_nigra.vcf")
```

----

### `PopGen.write`
```julia
PopGen.write(data::PopData; filename::String, kwargs...)
```
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the docstring to `genepop` with `?genepop`.

| File Format | Extensions             | Docstring            |
| :---------- | :--------------------- | :------------------- |
| genepop     | `.gen`, `.genepop`     | `?genepop`   |
| structure   | `.str`, `.structure`   | `?structure` |
| delimited   | `.csv`, `.txt`, `.tsv` | `?delimited` |

**Example**
```
cats = @nancycats;
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
