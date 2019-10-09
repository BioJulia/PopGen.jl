## Milligan 2002 dyadic relatedness - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1462494/pdf/12663552.pdf

cats=nancycats()

cats.samples

cat1=PopGen.get_genotype(cats, sample = "N111", locus = "fca78")
cat2=PopGen.get_genotype(cats, sample = "N100", locus = "fca78")

allele_freq = PopGen.allele_freq_mini(cats.loci.fca78)

function pr_l_s(x, y, allele_freq)
    ## Calculate Pr(Li | Sj)

    if x[1] == x[2] == y[1] == y[2] ##l1 - AiAi AiAi
        p = allele_freq[x[1]]
        [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]

    elseif (x[1] == x[2]) & (y[1] == y[2]) & (x[1] != y[1]) #L2 - AiAi AjAj
        p = (allele_freq[x[1]], allele_freq[y[1]])
        [0, prod(p), 0, prod(p) * p[2], 0, prod(p) * p[1], 0, 0, prod(p.^2)]

    elseif (x[1] == x[2] == y[1]) & (x[1] != y[2]) #L3 - AiAi AiAj
        p = (allele_freq[x[1]], allele_freq[y[2]])
        [0, 0, prod(p), 2 * prod(p) * p[1], 0, 0, 0, prod(p) * p[1], 2 * prod(p) * p[1]^2]

    elseif (x[1] == x[2]) & (y[1] != y[2]) & (x[1] != y[1]) & (x[1] != y[2]) #L4 - AiAi AjAk
        p = (allele_freq[x[1]], allele_freq[y[1]], allele_freq[y[2]])
        [0, 0, 0, 2 * prod(p), 0, 0, 0, 0, 2 * prod(p) * p[1]]

    elseif (x[1] == y[1] == y[2]) & (x[1] != x[2]) #L5 - AiAj AiAi
        p = (allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, prod(p), 2 * prod(p) * p[1], 0, prod(p) *p[1], 2 * prod(p) * p[1]^2]

    elseif (x[1] != x[2]) & (y[1] == y[2]) & (x[1] != y[1]) & (x[2] != y[1]) #L6 - AjAk AiAi
        p = (allele_freq[y[1]], allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, 0, 2 * prod(p), 0, 0, 2 * prod(p) * p[1]]

    elseif (x[1] == y[1]) & (x[2] == y[2]) & (x[1] != x[2]) #L7 - AiAj AiAj
        p = (allele_freq[x[1]], allele_freq[x[2]])
        [0, 0, 0, 0, 0, 0, 2 * prod(p), prod(p) * sum(p), 4 * prod(p.^2)]

    elseif (x[1] == y[1]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[2] != y[2]) #L8 - AiAj AiAk
        p = (allele_freq[x[1]], allele_freq[x[2]], allele_freq[y[2]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    elseif (x[1] != x[2]) & (x[1] != y[1]) & (x[1] != y[2]) & (x[2] != y[1]) & (x[2] != y[2]) & (y[1] != x[2]) #L9 - AiAj AkAl
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

dyad_likelihood(tst, rand(dirichlet_prior))
#Need to either maximize this sum or use it as the likelihood in a bayesian model and sample from the posterior.


#Need to condition on each of the 9 delta coefficients to maximize the likelihood and then calculate the relatedness

using Optim

lower = [0.0, 0, 0, 0, 0, 0, 0, 0, 0]
upper = [1.0, 1, 1, 1, 1, 1, 1, 1, 1]

dirichlet_init = Dirichlet(9, 1)

optimize(Δ -> dyad_likelihood(Δ, tst), lower, upper, rand(dirichlet_prior), Fminbox(NelderMead()))
## Not currently adhering to requirement that ΣΔ == 1




## Calculate theta and r
function relatedness_calc(Δ)
    θ = Δ[1] + 0.5 * (Δ[3] + Δ[5] + Δ[7]) + 0.25 * Δ[8]
    2 * θ
end
