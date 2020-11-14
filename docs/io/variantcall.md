---
id: vcf
title: Variant Call Format
sidebar_label: Variant Call Format
---

## Import a BCF/VCF file as `PopData`

```julia
vcf(infile::String; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
bcf(infile::String; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
```
PopGen.jl provides the commands `vcf` and `bcf` to import a variant call format files into `PopData`. The reader also accepts files that are gzipped. 

:::note Lazy Loading
The packages required to import BCF/VCF files (`GeneticVariation.jl` and `GZip.jl`) are lazy-loaded, therefore you must invoke `using GeneticVariation, GZip` to get this import functionality (see example below). 
:::


### Arguments

- `infile::String` : path to file, in quotes. **must end in `.gz` if gzipped**

### Keyword Arguments

- `rename_loci::Bool`: whether to simplify loci names to `snp_#` (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci (default: `false`)
- `silent::Bool` : whether to print file information during import (default: `false`)


### Example

```julia
using GeneticVariation, GZip
potato = vcf("/home/data/russet_potatoes.vcf.gz", allow_monomorphic = true)
cabbage = bcf("/home/data/nappa_cabbage.bcf", rename_loci = true, silent = true)
```

### Mixed-Ploidy data
In the event your variant call file is for mixed-ploidy data (where ploidy is not the same across all samples, e.g. PoolSeq), you will need to perform an additional step after reading in your data as `PopData` to convert the `.loci.genotype` column into a `GenoArray`:

```julia
julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;

julia> mydata.loci.genotype =  mydata.loci.genotype |> Array{Union{Missing, NTuple}}
```

### Format
Variant Call Format files follow a format standard, and while there is some wiggle-room for optional values, PopGen.jl only requires the core/mandatory components of a BCF/VCF, meaning problems should hopefully not arise regardless of which variant caller you are using (although we use `Freebayes` ourselves). Please open an issue if they do, or reach out to us on the community Slack.


### Allele encodings
When converting to `PopData`, the nucleotides will be recoded according to this table:

|    Base    |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |

:::caution Filter files beforehand
Keep in mind, BCF/VCF files need to be filtered **before** importing them into PopGen.jl. There is no and will be no VCF-filtering functionality to this package, as it is outside of the purpose of PopGen.jl. Refer to `vcftools`, `bcftools`, and `vcflib` to filter your sequence data. 
:::


### What BCF/VCF files contain and lack

Due to the nature of the file format, importing variant call files **will** provide:

- sample names
- ploidy of each sample
- locus names
- genotypes

but they **will not** provide:

- population information
- geographical coordinate information

This means you will need to add that information separately afterwards. Location data (which is *optional*) can be added to the `PopData` with the `locations!` command. Population names (*mandatory*) can be added using `populations!()`

## Acknowledgements
The majority of the BCF/VCF reader is thanks to the tremendous efforts of Ben J. Ward and the BioJulia contributors involved in [GeneticVariation.jl](https://github.com/BioJulia/GeneticVariation.jl), which we use to parse your files into `PopData` format. If you love the file importer, then give those folks your thanks. If something is wrong and/or you hate the importer, blame us first (and please [open up an issue](https://github.com/pdimens/PopGen.jl/issues) 😅). 