---
id: other_types
title: Other data types
sidebar_label: Other data types
---

While not strictly their own composite types, we also define aliases for genotypes and vectors of genotypes, as their explicit types can get a little unwieldy to use. The types shown below in the code blocks include their name and type (all types are of type `DataType`) on the first line, and what the alias actually defines on the second line.

### Genotype

```julia
Genotype::DataType
NTuple{N,<:Signed} where N
```

An `NTuple` is itself an alias for a `Tuple{Vararg{}}` , but you can think of it as Tuple of `N` length composed of items of a particular type, in this case it's items that are subtypes of `Signed` (the integer types). The length of the tuple (`N`) will vary based on the ploidy of the sample, and the element `Type` will vary whether the markers are snps (`Int8`) or microsatellites (`Int16`), making this a pretty flexible (but immutable) structure.

### GenoArray

```julia
GenoArray::DataType
AbstractVector{S} where S<:Union{Missing,Genotype}
```

As you can guess from the name, this Type wraps `Genotype` into a Vector, while permitting `missing` values (what's genetics without missing data!?). By using `AbstractVector` (rather than `Vector`), we also have the flexibility of functions working on things like `SubArrays` out of the box. 

:::note why bother defining these aliases?
Getting the most out of Julia and demonstrating good practices means making sure functions work on the things they're supposed to, and give informative error messages when the input isn't suitable for the function (a rare case of _wanting_ MethodErrors). Without these aliases, functions would either have vague definitions like `f(x,y,z) where x <: AbstractArray` and potentially cause errors, or overly complicated definitions like `f(x::AbstractVector{S},y,z) where {N, T<:Signed,S<:NTuple{N,T}}` and not be very legible. Instead, functions are written as `f(x,y,z) where x<:GenotypeArray`, and that seems like a good compromise of getting the latter while looking like the former.