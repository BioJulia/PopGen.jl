export read_from, file_import
#TODO update docs
"""
    read_from(infile::String; kwargs...)
Wraps `csv()`, `genepop()`, `bcf()`, and `vcf()` to read a file in as a `PopObj`. File type is
inferred from the file extension (case insensitive):
- delimited: .csv | .tsv | .txt
- genepop: .gen | .genepop
- variant call format: .vcf
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

    elseif ext in ["vcf", "bcf"]
        ext == "vcf" && return vcf(infile)
        ext == "bcf" && return bcf(infile)

    else
        @error "File type not recognized by filename extension \n delimited: .csv | .tsv | .txt \n genepop: .gen | .genepop \n variant call formant: .bcf | .vcf"
    end
end

const file_import = read_from

#TODO add to Read and API docs
"""
    write_to(data::PopData; filename::String, kwargs...)
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?popdata2filetype`. For example, to find the appropriate keywords for a converstion
to Genepop format, call up the docstring to `popdata2genepop` with `?popdata2genepop`.

### Example
```
cats = nancycats();
fewer_cats = omit_samples(cats, samples(cats)[1:10]);
write_to(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
### File format | extensions |  docstring
- genepop | `.gen` or `.genepop` | `?popdata2genepop`
"""
function write_to(data::PopData; filename::String, kwargs...)
    ext = split(filename, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        popdata2genepop(data, filename = filename, kwargs...)
    else
        @error "File type not recognized by filename extension. Please see the docstring"
    end
end