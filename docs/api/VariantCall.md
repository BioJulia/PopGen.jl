---
id: variantcall
title: VariantCall.jl
sidebar_label: VariantCall.jl
---
## VariantCall.jl

### `bcf`
```julia
bcf(infile::String)
```
Load a BCF file into memory as a PopData object. Population and [optional] location information need to be provided separately.
- `infile` : path to BCF file
- `allow_monomorphic` : whether to keep monomorphic loci in the dataset (default: `false`)
- `silent::Bool`   : whether to print file information during import (default: `true`)


### `vcf`
```julia
vcf(infile::String)
```
Load a VCF file into memory as a PoDataj object. Population and [optional] location information need to be provided separately.
- `infile` : path to VCF file
- `allow_monomorphic` : whether to keep monomorphic loci in the dataset (default: `false`)
- `silent::Bool`   : whether to print file information during import (default: `true`)
