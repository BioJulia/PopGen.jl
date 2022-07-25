---
id: readingdata
title: Read/Write data
sidebar_label: Read/Write data
---

PopGen.jl (via PopGenCore.jl) provides a handful of file readers and writers with which to create `PopData`. Each of the file types have their own file reader denoted simply by the file type:

| File Format         | Extensions             | Docstring     | Read | Write |
| :------------------ | :--------------------- | :------------ | :---:|:---:|
| [delimited](delimited.md)| `.csv`, `.txt`, `.tsv` | `?delimited`  |👍 | 👍 |
| [genepop](genepop.md)| `.gen`, `.genepop`     | `?genepop`    |👍 | 👍 |
| [structure](structure.md)| `.str`, `.structure`   | `?structure`  |👍 | 👍 |
| [plink](plink.md) (ped) | `.ped` | `?plink` | 👍 | 👍 |
| [variant call format (vcf)](variantcall.md) | `.vcf`, `.vcf.gz`| `?vcf`  |👍 | |
| [variant call format (bcf)](variantcall.md) | `.bcf`, `.bcf.gz`| `?bcf`  | 👍| |
| [baypass](baypass.md) | `.baypass` | `?baypass` | | 👍|


## Read in data
You're encouraged to use these functions, but PopGen.jl also provides you with an all-encompassing wrapper  `PopGen.read()`. Given the ubiquity of the function name, it is not exported. If using PopGenCore.jl directly, you will need to call it with `PopGenCore.read`.

:::caution Windows users
Make sure to change the backslashes `\` in your file path to double-backslashes `\\` or forward slashes `/` 
:::

:::note monomorphic loci
By default, the file reading methods drop monomorphic loci and inform you which were removed, so do not be alarmed if the number of loci in your `PopData` is different from the source data. You can disable this
behavior with the argument `allow_monomorphic = true`. Monomorphic loci are removed by default because they
can give spurious/misleading results for some analyses, such as kinship estimators.
:::

## `PopGen.read()`

```julia
PopGen.read(infile::String; kwargs...)
```

where `infile` is a String of your filename (in quotes) and `kwargs` are the corresponding keyword arguments associated with your file type. The function `PopGen.read()` uses all the same keyword arguments as do the commands specific to their file types, therefore you should have a look at those commands (usually the defaults suffice).

`PopGen.read()` infers the file type from the file extension, so for it to work properly your file must end with the extensions permitted below (case insensitive). If you're feeling particularly rebellious and your file does not conform to these extensions (such as a genepop file with a `.gen.final.v2.seriously` extension), then feel free to use the specific file importers, since they use the same exact syntax, there is zero difference in performance, and ignore file extensions. Ultimately, what crazy extensions you give your files is your business, and we love that about you.

**Examples**
```julia
salmon = PopGen.read("o_mykiss.gen", digits = 3, popsep = "SALMON")
ginko = PopGen.read("g_biloba.txt", delim = ",", digits = 2, silent = true)
```

## Write PopData to file

```julia
PopGen.write(data::PopData; filename::String, kwargs...)
```
To complement `PopGen.read()`, PopGen.jl offers `PopGen.write()`, 
which writes `PopData` to different file formats. Like the file 
reader, `PopGen.write()` will infer the correct output file type 
from the output filename's extensions. Given the ubiquity of the 
function name, it is not exported. If using PopGenCore.jl directly, 
you will need to call it with `PopGenCore.write`. 

Additional keyword arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific file writer with the format `?filetype` like shown above. For example, to find the appropriate keywords for a conversion to Genepop format, call up the docstring to `genepop` with `?genepop`.

**Examples**
```julia
cats = @nancycats;
fewer_cats = omit(cats, names = samplenames(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "horizontal")
PopGen.write(fewer_cats, filename = "filtered_nancycats.txt", digits = 4, format = "tidy", delim = ",")
```
