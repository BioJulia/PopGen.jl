---
slug: relatedness
title: Relatedness Tutorial
author: Pavel Dimens
author_title: Little this, little that
author_url: https://github.com/pdimens
author_image_url: https://avatars1.githubusercontent.com/u/19176506?s=460&u=3afad1d1ef3b09ddc4ab7108143f515be3412d5a&v=4
tags: [tutorials]
---

## Getting Started
In a population genetics study, you often need to identify if there are kin in your data. This may be necessary because you are trying to remove kin from your data (because of Hardy-Weinberg assumptions), or maybe kinship is a central interest in your study. Either way, the goal of this tutorial is to provide you with a basic tutorial on using PopGen.jl to perform a relatedness analysis, which is sometimes called a _kinship_ analysis. To follow along, you'll need to have Julia, along with the packages PopGen.jl, PopGenSims.jl, and StatsBase.jl installed. We'll be using the `nancycats` data because it's smaller than `gulfsharks`, so things should be a lot quicker.

```julia
using PopGen, PopGenSims, StatsBase

julia> cats = nancycats();
```

### Estimators
Like `Coancestry` and the R packages that wrap it (i.e. `relate`, `related`), PopGen.jl provides a whole bunch of relatedness estimators that you can choose from for your data. Unfortunately, there is no right answer and you will need to use your discretion. Some people choose an estimator based on the heterozygosity of thhe data, others choose one based on more liberal or conservative values, and there are yet more criteria one can consider for choosing an estimator. To keep things simple, we're going to use `LynchLi`. Why? Because I'm the one writing this tutorial, and I said so :grin:. 


## The Steps
### 1. Calculate pairwise relatedness
This seems pretty obvious, but it needs to be said. In order to do the analysis, you need you get the pairwise relatedness values for the individuals of interest in your data. To keep things simple, we're going to do an all-by-all comparison. But, we also want to boostrap the pairs to create confidence intervals ("CI") for each pair, so let's talk about that.

### 2. Bootstrap to calculate CI
It would behoove us to bootstrap the loci in a pairwise comparison _n_ number of times so we can create a confidence interval for the relatedness estimates for each pair. This inflates the runtime of the analysis substantially, but it's a very useful method in making sense of our estimates. If one was to be conservative (and we generally are), then we would reject an estimate for a pair whose CI includes zero. In PopGen.jl, the estimates and bootstrapping are done all at once to minimize processing time, so the command for that would be
```julia
julia> rel_out  = relatedness(cats, method = LynchLi, iterations = 1000)
```
By default, the `relatedness` function uses a 95% CI (`interval = (0.025, 0.975)`), but you can change that with `interval = (low,high)` where `low` and `high` are decimals of your quantiles. 

### 3. Create CI's for the sibships

:::note There's more???
The next set of steps seem like a lot more work, so allow me to explain. The estimators all generally give you some value between 0-1 (or 0-0.5, same idea) and you can intuit that certain values mean certain things, like that `0` is "unrelated", `0.25` is "half-sib", and `0.5` is "full-sib". However, those are fixed values, so how do we know how far we can deviate from 0.25 (for example) and still call our pair half-siblings? Instead of hand-waving, we can create confidence intervals from simulated data to act as sibship ranges for our data. If this doesn't make sense yet, it will below. Promise! 
:::

#### i. simulate known sibship pairs
Next, we need to further contextualize what out estimates actually mean, and we need to devise a rigorous and defensible method to define the sibling-ship ("sibship") of a pair of samples as **unrelated**, **half-sibs**, or **full-sibs**. To do this, we're going to use PopGenSims.jl to simulate data based on the allele frequencies in our data. What `simulate_sibship` does is simulate individuals based on the allele frequencies in your `PopData`, then simulate offspring of a particular siblingship resulting from the "mating" of those individuals. For example, when you specify `"fullsib"`, it generates two offspring resulting from two parents, `n` number of times. We want to do this for each of the three kinds of sibships.

```julia
julia> unrelated_sims = simulate_sibship(cats, n = 500, relationship= "unrelated")
julia> halfsib_sims = simulate_sibship(cats, n = 500, relationship= "halfsib")
julia> fullsib_sims = simulate_sibship(cats, n = 500, relationship= "fullsib")
```

Technically, we could merge all three results into a single `PopData`, but it will be easier to explain things if we don't.  

#### ii. get relatedness estimates for the simulated data
Next, we want to get the relatedness estimate for each simulated pair of "known" sibship. We are only interested in the values for the simulated pairs and not samples across pairs. If you aren't sure why that is, think of it this way: we're trying to create a range of values where we can confidently say unknown things are full-sibs (or half-sib, etc.), so we want to know what range of values we get from a bunch of known fullsib pairs, not the unknown relationships of samples between pairs. 

It would a nightmare to manually specify only 2 individuals at a time into `relatedness()`. Instead, the function has a shortcut built into it that will recognize the `population` names generated from `simulate_sibship` and only give you relatedness estimates for those pairs. So, we just need to run it once for each of our simulations, this time _without_ bootstrapping because we are only interested in the estimates. Make sure to use the same estimator! The run will be a lot faster this time because it only needs to calculate estimates for 500 pairs each time (our `n` from above) and without bootstrapping.

```julia
jula> un_sims_rel = relatedness(unrelated_sims, method = LynchLi)
jula> half_sims_rel = relatedness(halfsibs_sims, method = LynchLi)
jula> full_sims_rel = relatedness(fullsibs_sims, method = LynchLi)
```

Take a breath, we're almost there!

#### iii. sibship intervals
What we just did is create null distributions for each sibship relationship, so now all that's left is to get a confidence interval from each.

```julia
jula> unrelated_ci = quantile(un_sims_rel.LynchLi, (0.025, 0.975))
jula> halfsibs_ci = quantile(half_sims_rel.LynchLi, (0.025, 0.975))
jula> fullsibs_ci = quantile(full_sims_rel.LynchLi, (0.025, 0.975))
```

### 4. Final data assessment
Now that we have our relatedness estimates and the acceptable sibship ranges given our data, we can see where our data falls.
_Now_, we can say that a particular sample pair is unrelated/halfsib/fullsib if:
1. that pair's confidence interval does not include zero and 
2. that pair's estimate falls within any of the three calculate ranges


### Final remarks
There is more that can be done for relatedness, like network analysis. However, this tutorial covers what we consider a reasonable way to assess kinship in population genetic studies. If removing kin is your ultimate goal, consider the effects doing that may have on the analyses you are looking to do. Additionally, consider which individual of a pair you would remove and why. If you were curious as to how many identical loci a pair of samples may share, you can check that using `pairwise_identical()`. Good luck!