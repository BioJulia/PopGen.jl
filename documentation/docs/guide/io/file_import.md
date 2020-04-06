# Reading in data

Currently, PopGen.jl provides four different file parsers with which to create `PopData`. Each of the file types have their own file importer denoted simply by the file type:

| File type                     | Function to import it |
| ----------------------------- | --------------------- |
| Delimited                     | `delimited()`         |
| Genepop                       | `genepop()`           |
| Variant Call Format           | `vcf()`               |
| Variant Call Format (bzipped) | `bcf()`               |

You're encouraged to use functions, but PopGen.jl also provides you with an all-encompassing wrapper  called `read_in()`. 

::: details file_import alias
The name "read_in" is short and sweet, but isn't the most descriptive, so we provide a more explicitly-named alias `file_import()`. Feel free to use whichever you like best.
:::

```julia
read_in(infile::String; kwargs...)
file_import(infile::String; kwargs...)
```

where `infile` is a String of your filename (in quotes) and `kwargs` are the corresponding keyword arguments associated with your file type. The function `read_in()` uses all the same keyword arguments as do the commands specific to their file types, therefore you should have a look at those commands (usually the defaults suffice). 

`read_in()` infers the file type from the file extension, so for it to work properly your file must end with the extensions permitted below (case insensitive). If you're feeling particularly rebellious and your file does not conform to these extensions (such as a genepop file with a `.gen.final.v2.seriously` extension), then feel free to use the specific file importers, since they use the same exact syntax, there is zero difference in performance, and ignore file extensions. Ultimately, what crazy extensions you give your files is your business, and we love that about you. 

## Supported File Types

### [Delimited files](delimited.md) 

Accepted extensions: `.csv`, `.tsv`, `.txt`

- files in which values are separated using a consistent delimiter, such as commas, spaces, or tabs
- first rows are column names
- each line represents a row


### [Genepop Files](genepop.md)

Accepted extensions: `.gen`, `.genepop`

- first row is a comment and skipped
- then comes list of all loci, usually 1-per-line
  - sometimes horizontally arranged and separated by commas
- populations separated by a word like `"POP"`
- sample names followed by a **comma, then a tab or space**
- genotypes separated by **tabs** or **spaces**
- genotypes represented as a combination of ploidy x _n_-digits 
	- e.g. for genotype 001002 (3 digits per allele)
	- allele 1 = 001
	- allele 2 = 002


### [Variant Call Format](variantcall.md)

Accepted extensions: `.vcf`, `.bcf`

This format is **much** more complex and variable depending on which variant caller produced the file. If you're super duper interested in the specifications of BCF/VCF files, have a look at the [official specification documentation](http://samtools.github.io/hts-specs/VCFv4.3.pdf).