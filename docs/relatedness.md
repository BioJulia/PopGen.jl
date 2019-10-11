## Background

Sometimes you want or need to know the relatedness of individuals in a sample. Relatedness is exactly what its name implies: how related individuals are given some provided genetic information. Relatedness can be used in quantitative genetics to estimate heritability, additive genetic variances, and covariances. It can also be used in population genetics to study isolation-by-distance or population structure. 

Given two diploid individuals, there are 9 "identity by descent" models available between them ([Jacquard 1975](https://www.springer.com/gp/book/9783642884177), paywall), as shown below (from [Milligan 2003](https://www.genetics.org/content/163/3/1153.full)): 

[![Jacquard IBD](img/jacquard_identitiies.jpg)](https://www.genetics.org/content/163/3/1153.full) 

This means that for large not-inbred populations, we can reduce this to the probability a pair of individuals share two alleles (by descent), and the probability they share one. There have been a handful of brilliant minds who have come up with methods for estimating relatedness via these relationships, and currently PopGen.jl implements the maximum likelihood approach presented in [Milligan 2003](https://www.genetics.org/content/163/3/1153.full).  You can imagine there's a lot that  happens under the hood to perform this for all loci across all individuals-- all of which dutifully written by :star2: Jason Selwyn :star2:. 



## Calculate Relatedness

The docstring for the function will go here. This is filler text for now. Mentioning of taking ___ and ____ and input and returns a _______. 

```julia
relatedness(data::PopObj, arg::Type, arg::Type, arg::Type)
```

### arguments

- 
- 
- 

### keyword arguments

- 
- 


### example

```julia tab="relatedness"

```

```tab="output"

```

See [Development](hidden_api.md) for all of the APIs associated with `relatedness()` 

## Credits

Both [Convex.jl](https://github.com/JuliaOpt/Convex.jl) and [ECOS.jl](https://github.com/JuliaOpt/ECOS.jl) are pivotal for these calculations, and we thank the authors for their time developing and maintaining them, along with the members of the Julia Slack channel for pointing us towards those packages.