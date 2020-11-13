export read_from, file_import, write_to

#TODO update docs with edits to this docstring and structure info
"""
    read_from(infile::String; kwargs...)
Wraps the individual file importers to read a file in as a `PopData` object. File type is
inferred from the file extension (case insensitive): \n

| File Format         | Extensions             | Docstring     |
| :------------------ | :--------------------- | :------------ |
| delimited           | `.csv`, `.txt`, `.tsv` | `?delimited`  |
| genepop             | `.gen`, `.genepop`     | `?genepop`    |
| structure           | `.str`, `.structure`   | `?structure`  |
| variant call format (vcf) | `.vcf`, `.vcf.gz`| `?vcf`  |
| variant call format (bcf) | `.bcf`, `.bcf.gz`| `?bcf`  |

This function uses the same keyword arguments (and defaults) as the file importing
functions it wraps; please see their respective docstrings in the Julia help console.
(e.g. `?genepop`) for specific usage details. Use the alias function `file_import`
interchangeably if you prefer the explicit name instead.

## Examples
```
read_from("cavernous_assfish.gen", digits = 3)

file_import("bos_tauros.csv", silent = true)

read_from("juglans_nigra.vcf")
```
"""
function read_from(infile::String; kwargs...)
    ext = split(infile, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        return genepop(infile;kwargs...)

    elseif ext in ["csv", "txt", "tsv"]
        return delimited(infile; kwargs...)

    elseif ext in ["str", "structure"]
        return structure(infile; kwargs...)

    elseif ext in ["vcf", "bcf"]
        ext == "vcf" && return vcf(infile)
        ext == "bcf" && return bcf(infile)

    else
        @error "File type not recognized by filename extension \n delimited: .csv | .tsv | .txt \n genepop: .gen | .genepop \n structure: .str | .structure \n variant call formant: .bcf | .vcf"
    end
end

const file_import = read_from

"""
    write_to(data::PopData; filename::String, kwargs...)
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?popdata2filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the docstring to `popdata2genepop` with `?popdata2genepop`.

| File Format | Extensions             | Docstring          |
| :---------- | :--------------------- | :----------------- |
| genepop     | `.gen`, `.genepop`     | ?popdata2genepop   |
| JLD2        | `.jld2`                | ?popdata2jld2      |
| delimited   | `.csv`, `.txt`, `.tsv` | ?popdata2delimited |
| structure   | `.str`, `.structure`   | ?popdata2structure |

### Example
```
cats = nancycats();
fewer_cats = omit(cats, name = samples(cats)[1:10]);
write_to(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
"""
function write_to(data::PopData; filename::String, kwargs...)
    ext = split(filename, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        popdata2genepop(data, filename = filename; kwargs...)
    elseif ext == "jld2"
        popdata2jdl2(data, filename = filename)
    elseif ext in ["str", "structure"]
        return popdata2structure(data; kwargs...)
    elseif ext in ["csv", "txt", "tsv"]
        popdata2delimited(data, filename = filename; kwargs...)
    else
        @error "File type not recognized by filename extension. Please see the docstring"
    end
end


"""
    popdata2jdl2(data::PopData; filename::String)
Write PopData to a `JLD2` file.
"""
function popdata2jdl2(data::PopData; filename::String)
    JLD2.@save filename data
end
