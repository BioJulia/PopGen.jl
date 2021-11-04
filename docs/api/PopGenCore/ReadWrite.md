---
id: read
title: ReadWrite.jl
sidebar_label: ReadWrite.jl
---
## PopGenCore.jl/src/io/ReadWrite.jl
📦  => not exported | 
🟪 => exported by PopGenCore.jl | 
🔵 => exported by PopGen.jl

### 📦 read
```julia
PopGen.read(infile::String; kwargs...)
```
Wraps the individual file importers to read a file in as a `PopData` object. File type is
inferred from the file extension (case insensitive):

| File Format         | Extensions             | Docstring     |
| :------------------ | :--------------------- | :------------ |
| delimited           | `.csv`, `.txt`, `.tsv` | `?delimited`  |
| genepop             | `.gen`, `.genepop`     | `?genepop`    |
| structure           | `.str`, `.structure`   | `?structure`  |
| plink               | `.bed`, `.ped`  | `?plink`  |
| variant call format (vcf) | `.vcf`, `.vcf.gz`| `?vcf`  |
| variant call format (bcf) | `.bcf`, `.bcf.gz`| `?bcf`  |

This function uses the same keyword arguments (and defaults) as the file importing
functions it wraps; please see their respective docstrings in the Julia help console.
for specific usage details (e.g. `?genepop`).

**Examples**
```
PopGen.read("cavernous_assfish.gen", digits = 3)
PopGen.read("juglans_nigra.vcf")
```

----

### 📦 write
```julia
PopGen.write(data::PopData, filename::String, kwargs...)
PopGen.write(data::PopData; filename::String, kwargs...)
```
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the `?genepop` docstring.

| File Format | Extensions             | Docstring          |
| :---------- | :--------------------- | :----------------- |
| genepop     | `.gen`, `.genepop`     | ?genepop   |
| delimited   | `.csv`, `.txt`, `.tsv` | ?delimited |
| structure   | `.str`, `.structure`   | ?structure |

**Example**
```
cats = @nancycats;
fewer_cats = omit(cats, name = samplenames(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```