---
id: variantcall
title: VariantCall.jl
sidebar_label: VariantCall.jl
---

### `openvcf`
    openvcf(::String)
Open VCF file (`.vcf/.gz`, or `.bcf/.gz`) and return an `IO` stream in reading mode `"r"`.
Adapted from OpenMendel/VCFTools.jl
https://github.com/OpenMendel/VCFTools.jl/blob/master/src/gtstats.jl#L169


----

### `bcf`
    bcf(infile::String; ; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
Load a BCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to BCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:

|Base| | Allele |
|:---:|:---:|:---:|
| A | => | 1 |
| T | => | 2 |
| C | => | 3 |
| G | => | 4 |

#### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;
julia> mydata.loci.genotype =  mydata.loci.genotype |> Array{Union{Missing, NTuple}}
```
----

### `vcf`
    vcf(infile::String; ; rename_snp::Bool, silent::Bool, allow_monomorphic::Bool)
Load a VCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to VCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:

|Base| | Allele |
|:---:|:---:|:---:|
| A | => | 1 |
| T | => | 2 |
| C | => | 3 |
| G | => | 4 |
#### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = vcf("path/to/file.vcf", silent = true, rename_loci = true) ;
julia> mydata.loci.genotype =  mydata.loci.genotype |> Array{Union{Missing, NTuple}}
```