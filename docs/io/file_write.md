---
id: file_write
title: Writing data to file
sidebar_label: Writing data
---

```julia
PopGen.write(data::PopData; filename::String, kwargs...)
```
To complement `PopGen.read()`, PopGen.jl offers `PopGen.write()`, 
which writes `PopData` to different file formats. Like the file 
reader, `PopGen.write()` will infer the correct output file type 
from the output filename's extensions. Given the ubiquity of the 
function name, it is not exported. If using PopGenCore.jl directly, 
you will need to call it with `PopGenCore.write`. 
PopGen.jl supports writing to these file formats:

| File Format | Extensions             | Docstring            |
| :---------- | :--------------------- | :------------------- |
| genepop     | `.gen`, `.genepop`     | `?genepop`   |
| delimited   | `.csv`, `.txt`, `.tsv` | `?delimited` |
| Structure/fastStructure   | `.str`, `.structure` | `?structure` |

Additional keyword arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific file writer with the format `?filetype` like shown above. For example, to find the appropriate keywords for a conversion to Genepop format, call up the docstring to `genepop` with `?genepop`.

** Examples **
```julia
cats = @nancycats;
fewer_cats = omit(cats, names = samples(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "horizontal")
PopGen.write(fewer_cats, filename = "filtered_nancycats.txt", digits = 4, format = "tidy", delim = ",")
```
