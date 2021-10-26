#Run all lines of code in PopGen.jl
#Then read in all functions within PairwiseRelatedness.jl

verbose = true
data = @nancycats

alleles = Dict()
for loc in loci(data)
    alleles[loc] = allelefreq(locus(data, loc))
end
alleles

#Working pair
ind1 = "N100"
ind2 = "N106"

#Not working pair ex.
ind1 = "N111"
ind2 = "N83"

L = all_loci_Pr_L_S(data, ind1, ind2, alleles)


Δ = Variable(8)
problem = maximize(sum(log(L * vcat(1 - sum(Δ), Δ))))
problem.constraints += 0 <= 1 - sum(Δ)
problem.constraints += 1 - sum(Δ) <= 1
problem.constraints += 0 <= Δ[1:8]
problem.constraints += Δ[1:8] <= 1

#shifted from actual relatedness calculations because the 1 - sum(Δ) goes at beginning
problem.constraints += 2 * ((1 - sum(Δ)) + 0.5 * (Δ[2] + Δ[4] + Δ[6]) + 0.25 * Δ[7]) <= 1
problem.constraints += 0 <= 2 * ((1 - sum(Δ)) + 0.5 * (Δ[2] + Δ[4] + Δ[6]) + 0.25 * Δ[7])

Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit = 100), verbose = verbose)



n = Variable(9)
problem2 = minimize(sum(log(L * exp(n)) / sum(exp(n))))
Convex.solve!(problem2, ECOSSolver())
