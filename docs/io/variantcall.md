---
id: vcf
title: Variant Call Format
sidebar_label: Variant Call Format
---

## Import a BCF/VCF file as `PopData`

```julia
vcf(infile::String; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
bcf(infile::String; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
```
PopGen.jl provides the commands `vcf` and `bcf` to import a variant call format files into `PopData`. The reader also accepts files that are gzipped. 

### Arguments

- `infile::String` : path to file, in quotes. **must end in `.gz` if gzipped**

### Keyword Arguments

- `rename_loci::Bool`: whether to simplify loci names to `snp_#` (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci (default: `false`)
- `silent::Bool` : whether to print file information during import (default: `false`)


### Example

```julia
cabbage = bcf("/home/data/nappa_cabbage.bcf", rename_loci = true, silent = true)
potato = vcf("/home/data/russet_potatoes.vcf.gz", allow_monomorphic = true)
```

### Mixed-Ploidy data
In the event your variant call file is for mixed-ploidy data (where ploidy is not the same across all samples, e.g. PoolSeq), you will need to perform an additional step after reading in your data as `PopData` to convert the `.genodata.genotype` column into a `GenoArray`:
```julia
julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;

julia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}
```

:::caution WIP
The extra step required by mixed-ploidy data is a work in progress. Feel free to submit a PR if you have ideas!
:::


### Format
Variant Call Format files follow a format standard, and while there is some wiggle-room for optional values, PopGen.jl only requires the core/mandatory components of a BCF/VCF, meaning problems should hopefully not arise regardless of which variant caller you are using (although we use `Freebayes` ourselves). Please open an issue if they do, or reach out to us on the community Slack.


### Allele encodings
When converting to `PopData`, the nucleotides will be recoded according to the table below. Note that this system differs slightly from
how PGDSpider2 recodes alleles (the 3 and 4 are switched).

|    Base    |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |

<details>
<summary>VCF Filtering</summary>

Keep in mind, BCF/VCF files need to be filtered **before** importing them into PopGen.jl. There is no and will be no VCF-filtering functionality to this package, as it is outside of the purpose of PopGen.jl. Refer to `vcftools`, `bcftools`, and `vcflib` to filter your sequence data. 

</details>

### What BCF/VCF files contain

Due to the nature of the file format, importing variant call files **will** provide:

- sample names
- ploidy of each sample
- locus names
- genotypes

### What BCF/VCF files lack
- population information
- geographical coordinate information

This means you will need to add that information separately afterwards. Location data (which is *optional*) can be added to the `PopData` with the `locations!` command. Population names (*mandatory*) can be added using `populations!()`

## Acknowledgements
The heavy lifting of the BCF/VCF reader is thanks to the tremendous efforts of the contributors involved with 
[GeneticVariation.jl](https://github.com/BioJulia/GeneticVariation.jl), and its successor [VariantCallFormat.jl](https://github.com/rasmushenningsson/VariantCallFormat.jl) 
which we use to parse files into `PopData` format. More specifically, the two packages use a file parser created from [Automa.jl](https://github.com/BioJulia/Automa.jl). If you love the file importer, then give those folks your thanks. If something is wrong and/or you hate the importer, blame us 
first (and please [open up an issue](https://github.com/biojulia/PopGenCore.jl/issues) 😅).
