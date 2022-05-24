---
id: fstglobal
title: FstGlobal.jl
sidebar_label: FstGlobal.jl
---
## PopGen.jl/src/FStatistics/FstGlobal.jl
| 📦  not exported | 🔵  exported by PopGen.jl |
|:---:|:---:|

### 🔵 Hudson
```julia
Hudson(data::PopData)
Hudson(pop1::T, pop2::T) where T<:GenoArray
Hudson(population_1::T, population_2::T) where T<:AbstractMatrix
```

----

### 🔵 Nei
```julia
Nei(data::PopData)
Nei(population_1::T, population_2::T) where T<:AbstractMatrix
```

----

### 🔵 WeirCockerham
```julia
WeirCockerham(data::PopData)
WeirCockerham(population_1::T, population_2::T) where T<:AbstractMatrix
```