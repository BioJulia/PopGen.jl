---
id: file_write
title: Writing data to file
sidebar_label: Writing data
---

```julia
write_to(data::PopData; filename::String, kwargs...)
```
To complement `read_from()`, PopGen.jl offers `write_to()`, which writes `PopData` to different file formats. Like the file reader, `write_to()`
will infer the correct output file type from the output filename's extensions. Currently, PopGen.jl supports writing to these file formats:

| File Format | Extensions             | Docstring            |
| :---------- | :--------------------- | :------------------- |
| genepop     | `.gen`, `.genepop`     | `?genepop`   |
| delimited   | `.csv`, `.txt`, `.tsv` | `?delimited` |
| Structure/fastStructure   | `.str`, `.structure` | `?structure` |

Additional keyword arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific file writer with the format `?filetype` like shown above. For example, to find the appropriate keywords for a conversion to Genepop format, call up the docstring to `genepop` with `?genepop`.

```julia
cats = @nancycats;
fewer_cats = omit(cats, names = samples(cats)[1:10]);
write_to(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "horizontal")
write_to(fewer_cats, filename = "filtered_nancycats.txt", digits = 4, format = "tidy", delim = ",")
```
