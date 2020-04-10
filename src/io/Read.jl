export read_in, file_import

"""
    read_in(infile::String; kwargs...)
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
read_in("cavernous_assfish.gen", digits = 3)

file_import("bos_tauros.csv", silent = true)

read_in("juglans_nigra.vcf")
```
"""
function read_in(infile::String; kwargs...)
    ext = split(infile, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        return genepop(infile;kwargs...)

    elseif ext in ["csv", "txt", "tsv"]
        return delimited(infile; kwargs...)

    elseif ext in ["vcf", "bcf"]
        ext == "vcf" && return vcf(infile)
        ext == "bcf" && return bcf(infile)

    else
        @error "file type not recognized by filename extension \n delimited: .csv | .tsv | .txt \n genepop: .gen | .genepop \n variant call formant: .bcf | .vcf"
    end
end

const file_import = read_in
