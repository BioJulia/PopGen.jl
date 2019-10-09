## Milligan 2002 dyadic relatedness - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1462494/pdf/12663552.pdf

cats=nancycats()

cats.samples

cat1=PopGen.get_genotype(cats, sample = "N128", locus = "fca23")
cat2=PopGen.get_genotype(cats, sample = "N100", locus = "fca23")

allele_freq = PopGen.allele_freq_mini(cats.loci.fca23)

function pr_l_s(x, y, allele_freq)
    #= Current Bugs
    1) For the genotypes 138|146 & 136|146 it is not properly sorted into the correct class (L8) because of the order of the alleles

    =#

    #= Calculate Pr(Li | Sj)
    If the allele identity falls into this class (L1-L9), generate the
    probabilities of it belonging to each of the different classes and
    return that array of 9 distinct probabilities
    =#

    ## class L1 -  AᵢAᵢ AᵢAᵢ ##
    if x[1] == x[2] == y[1] == y[2] ##l1 -
        p = allele_freq[x[1]]
        [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]

    ## class L2 - AᵢAᵢ AⱼAⱼ ##
    elseif (x[1] == x[2]) & (y[1] == y[2]) & (x[1] != y[1])
        p = (allele_freq[x[1]], allele_freq[y[1]])
        [0, prod(p), 0, prod(p) * p[2], 0, prod(p) * p[1], 0, 0, prod(p.^2)]

    ## class L3 - AᵢAᵢ AᵢAⱼ ##
    elseif (x[1] == x[2] == y[1]) & (x[1] != y[2])
        p = (allele_freq[x[1]], allele_freq[y[2]])
        [0, 0, prod(p), 2 * prod(p) * p[1], 0, 0, 0, prod(p) * p[1], 2 * prod(p) * p[1]^2]

    ## class L4 - AᵢAᵢ AⱼAₖ ##
    elseif (x[1] == x[2]) & (y[1] != y[2]) & (x[1] != y[1]) & (x[1] != y[2])
        p = (allele_freq[x[1]], allele_freq[y[1]], allele_freq[y[2]])
        [0, 0, 0, 2 * prod(p), 0, 0, 0, 0, 2 * prod(p) * p[1]]

    ## L5 - AiAj AiAi ##
    elseif (x[1] == y[1] == y[2]) & (x[1] != x[2])
        p = (allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, prod(p), 2 * prod(p) * p[1], 0, prod(p) *p[1], 2 * prod(p) * p[1]^2]

    ## L6 - AjAk AiAi ##
    elseif (x[1] != x[2]) & (y[1] == y[2]) & (x[1] != y[1]) & (x[2] != y[1])
        p = (allele_freq[y[1]], allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, 0, 2 * prod(p), 0, 0, 2 * prod(p) * p[1]]

    ## L7 - AiAj AiAj ##
    elseif (x[1] == y[1]) & (x[2] == y[2]) & (x[1] != x[2])
        p = (allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, 0, 0, 2 * prod(p), prod(p) * sum(p), 4 * prod(p.^2)]

    ## L8 - AiAj AiAk ##
    elseif (x[1] == y[1]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[2] != y[2])
        p = (allele_freq[x[1]], allele_freq[x[2]], allele_freq[y[2]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    ## L9 - AiAj AkAl ##
    elseif (x[1] != x[2]) & (x[1] != y[1]) & (x[1] != y[2]) & (x[2] != y[1]) & (x[2] != y[2]) & (y[1] != x[2])
        p = (allele_freq[x[1]], allele_freq[x[2]], allele_freq[y[1]], allele_freq[y[2]])
        [0, 0, 0, 0, 0, 0, 0, 0, 4 * prod(p)]
    else
        [-9, -9, -9, -9, -9, -9, -9, -9, -9]
    end
end
tst = pr_l_s(cat1, cat2, allele_freq)

## Calculate Δ coefficients

function dyad_likelihood(Δ::Vector{Float64}, Pr_L_S::Vector{Float64})
    #Δ is what needs to be optimized
    #consist of 9 values between 0 and 1 which must also sum to 1
    #is then used to calculate relatedness

    -1 * sum(Pr_L_S .* Δ)
end


lower = [0.0, 0, 0, 0, 0, 0, 0, 0, 0]
upper = [1.0, 1, 1, 1, 1, 1, 1, 1, 1]

dirichlet_distr = Dirichlet(9, 1)
dyad_likelihood(rand(dirichlet_distr), tst)

optimized_Δ = optimize(Δ -> dyad_likelihood(Δ, tst), lower, upper, rand(dirichlet_distr), Fminbox(NelderMead()))

using JuMP
using GLPK
##### 8 and the 9th is 1-sum - here it is
function Δ_likelihood(Pr_L_S::Vector{Float64})
    #Δ is what needs to be optimized
    #consist of 9 values between 0 and 1 which must also sum to 1
    #is then used to calculate relatedness
    #push!(Δ, 1 - sum(Δ))

    model = Model(with_optimizer(GLPK.Optimizer))
    @variable(model, 0 <= Δ[1:8] <= 1)
    @objective(model, Max, sum(tst[1:8] .* Δ) + ((1 - sum(Δ)) * tst[9]))
    @constraint(model, con, 0 <= (1 - sum(Δ)) <= 1)
    optimize!(model)

    out = value.(Δ)
    push!(out, 1 - sum(out))

end

Δ_likelihood(tst)

#Need to either maximize this sum or use it as the likelihood in a bayesian model and sample from the posterior.


#Need to condition on each of the 9 delta coefficients to maximize the likelihood and then calculate the relatedness

lower = [0.0, 0, 0, 0, 0, 0, 0, 0, 0]
upper = [1.0, 1, 1, 1, 1, 1, 1, 1, 1]

dirichlet_distr = Dirichlet(9, 1)
dyad_likelihood(rand(dirichlet_distr), tst)


Δ[7] = 1.00
Δ[8] = 1.00


optimized_Δ = optimize(Δ -> Δ_likelihood(Δ, tst), lower, upper, rand(dirichlet_distr), Fminbox(NelderMead()))
## Not currently adhering to requirement that ΣΔ == 1
optimized_Δ.minimizer
optimized_Δ.iterations
optimized_Δ.initial_x


## Calculate theta and r
function relatedness_calc(Δ)
    θ = Δ[1] + 0.5 * (Δ[3] + Δ[5] + Δ[7]) + 0.25 * Δ[8]
    2 * θ
end

relatedness_calc(optimized_Δ.minimizer/sum(optimized_Δ.minimizer))
