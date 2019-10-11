#= Steps to finish:

    Solve dyadic optimization issues - set max higher?
        combined backtracking failure - need to isolate and find cause
            #N111 N83
            #N63 N64
        Relax tolerance?
    Solve dyadic inbreeding issues
    Output Δ coefficents

=#

#= Would be good to include

    Implement alternative relatedness metrics

=#

"""
    pr_l_s(x::Tuple, y::Tuple, alleles::Dict)
Calculate the probability of observing the particular allele state given each of
the 9 Jacquard Identity States for a single locus to create Table 1 from
Milligan 2002.
"""
function pr_l_s(x::Tuple, y::Tuple, alleles::Dict)
    #= Example
    cats = nancycats()
    cat1=get_genotype(cats, sample = "N100", locus = "fca23")
    cat2=get_genotype(cats, sample = "N111", locus = "fca23")
    allele = allele_freq(cats.loci.fca23)
    pr_l_s(cat1, cat2, allele)
    =#
    #=
    Calculate Pr(Li | Sj)
    If the allele identity falls into this class (L1-L9), generate the
    probabilities of it belonging to each of the different classes and
    return that array of 9 distinct probabilities
    =#

    ## class L1 -  AᵢAᵢ AᵢAᵢ ##
    if x[1] == x[2] == y[1] == y[2]
        p = alleles[x[1]]
        [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]

    ## class L2 - AᵢAᵢ AⱼAⱼ ##
    elseif (x[1] == x[2]) & (y[1] == y[2]) & (x[1] != y[1])
        p = (alleles[x[1]], alleles[y[1]])
        [0, prod(p), 0, prod(p) * p[2], 0, prod(p) * p[1], 0, 0, prod(p) * prod(p)]

    ## class L3a - AᵢAᵢ AᵢAⱼ ## - has issues because of allele order
    elseif ((x[1] == x[2] == y[1]) & (x[1] != y[2]))
        p = (alleles[x[1]], alleles[y[2]])
        [0, 0, prod(p), 2 * prod(p) * p[1], 0, 0, 0, prod(p) * p[1], 2 * prod(p) * p[1]^2]

    ## class L3b - AᵢAᵢ AⱼAᵢ ## - has issues because of allele order
    elseif ((x[1] == x[2] == y[2]) & (x[1] != y[1]))
        p = (alleles[x[1]], alleles[y[1]])
        [0, 0, prod(p), 2 * prod(p) * p[1], 0, 0, 0, prod(p) * p[1], 2 * prod(p) * p[1]^2]

    ## class L4 - AᵢAᵢ AⱼAₖ ##
    elseif (x[1] == x[2]) & (y[1] != y[2]) & (x[1] != y[1]) & (x[1] != y[2])
        p = (alleles[x[1]], alleles[y[1]], alleles[y[2]])
        [0, 0, 0, 2 * prod(p), 0, 0, 0, 0, 2 * prod(p) * p[1]]

    ## L5a - AiAj AiAi ## - has issues because of allele order
    elseif ((x[1] == y[1] == y[2]) & (x[1] != x[2]))
        p = (alleles[x[1]], alleles[x[2]])
        [0, 0, 0, 0, prod(p), 2 * prod(p) * p[1], 0, prod(p) *p[1], 2 * prod(p) * p[1]^2]

    ## L5b - AjAi AiAi ## - has issues because of allele order
    elseif (x[2] == y[1] == y[2] * (x[1] != x[2]))
        p = (alleles[x[2]], alleles[x[1]])
        [0, 0, 0, 0, prod(p), 2 * prod(p) * p[1], 0, prod(p) *p[1], 2 * prod(p) * p[1]^2]

    ## L6 - AjAk AiAi ##
    elseif (x[1] != x[2]) & (y[1] == y[2]) & (x[1] != y[1]) & (x[2] != y[1])
        p = (alleles[y[1]], alleles[x[1]], alleles[x[2]])
        [0, 0, 0, 0, 0, 2 * prod(p), 0, 0, 2 * prod(p) * p[1]]

    ## L7 - AiAj AiAj ##
    elseif (x[1] == y[1]) & (x[2] == y[2]) & (x[1] != x[2])
        p = (alleles[x[1]], alleles[x[2]])
        [0, 0, 0, 0, 0, 0, 2 * prod(p), prod(p) * sum(p), 4 * prod(p) * prod(p)]

    ## L8a - AiAj AiAk ##  - has issues because of allele order
    elseif ((x[1] == y[1]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[2] != y[2]))
        p = (alleles[x[1]], alleles[x[2]], alleles[y[2]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    ## L8b - AjAi AkAi ##  - has issues because of allele order
    elseif ((x[2] == y[2]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[1] != y[1]))
        p = (alleles[x[2]], alleles[x[1]], alleles[y[1]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    ## L8c - AjAi AiAk ##  - has issues because of allele order
    elseif ((x[2] == y[1]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[1] != y[2]))
        p = (alleles[x[2]], alleles[x[1]], alleles[y[2]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    ## L8d - AiAj AkAi ##  - has issues because of allele order
    elseif ((x[1] == y[2]) & (x[1] != x[2]) & (y[1] != y[2]) & (x[1] != y[1]))
        p = (alleles[x[1]], alleles[x[2]], alleles[y[1]])
        [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]

    ## L9 - AiAj AkAl ##
    elseif (x[1] != x[2]) & (x[1] != y[1]) & (x[1] != y[2]) & (x[2] != y[1]) & (x[2] != y[2]) & (y[1] != x[2])
        p = (alleles[x[1]], alleles[x[2]], alleles[y[1]], alleles[y[2]])
        [0, 0, 0, 0, 0, 0, 0, 0, 4 * prod(p)]
    else
        [-9, -9, -9, -9, -9, -9, -9, -9, -9]
    end
end


"""
    all_loci_Pr_L_S(data::PopObj, ind1::String, ind2::String, alleles::Dict)
Calculate the probability of observing the particular allele state given each of
the 9 Jacquard Identity States for all loci Creates Table 1 from Milligan 2002
"""
function all_loci_Pr_L_S(data::PopObj, ind1::String, ind2::String, alleles::Dict)
    #Need to skip loci where one or both individuals have missing data
    Pr_L_S = []
    for locus in String.(names(data.loci))
        #Extract the pair of interest's genotypes
        gen1 = get_genotype(data, sample = ind1, locus = locus)
        gen2 = get_genotype(data, sample = ind2, locus = locus)

        if gen1 !== missing && gen2 !== missing
            tmp = pr_l_s(gen1, gen2, alleles[locus])
            push!(Pr_L_S, tmp)
        end
    end
    return transpose(hcat(Pr_L_S...))
end

#Pr_L_S_inbreeding = all_loci_Pr_L_S(ind1, ind2, data, allele_frequencies)


"""
    no_inbreeding(Pr_L_S::LinearAlgebra.Transpose{Float64,Array{Float64,2}})
Remove Jacquard States which can only be the result of inbreeding
"""
function no_inbreeding(Pr_L_S::Transpose{Float64,Array{Float64,2}})
    for i in 1:6
        Pr_L_S[:,i] = 0 .* Pr_L_S[:,i]
    end
    return Pr_L_S
end

#Pr_L_S_noinbreeding = dyadic_ML(data, allele_frequencies) |> no_inbreeding


## Calculate Δ coefficients
#Need to either maximize this sum or use it as the likelihood in a bayesian model and sample from the posterior.
#currently only maximum likelihood optimization


"""
    Δ_optim(Pr_L_S::Transpose{Float64,Array{Float64,2}}, verbose::Bool)
Takes the probability of the allelic state given the identity by descent from
all available loci (either allowing for inbreeding or not) and calculated the
maximum likelihood Δ coefficients
"""
function Δ_optim(Pr_L_S::Transpose{Float64,Array{Float64,2}}, verbose::Bool = true)
    #Δ is what needs to be optimized
    #consist of 9 values between 0 and 1 which must also sum to 1
    #is then used to calculate relatedness

    Δ = Variable(9)
    problem = maximize(sum(log(Pr_L_S * Δ)))
    problem.constraints += Δ[9] == 1 - sum(Δ[1:8])
    problem.constraints += sum(Δ) == 1
    problem.constraints += 0 <= Δ[1:9]
    problem.constraints += Δ[1:9] <= 1
    #Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit=100, feastol=1e-7), verbose = verbose)
    #Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit=100, feastol=5e-6), verbose = verbose)
    Convex.solve!(problem, SCSSolver(verbose = verbose), verbose = verbose)

    Δ.value, problem.status
    #Should probably include some output that confirms that it did in fact converge and/or use multiple random starts to confirm not a local maxima
end

dyadicML_relatedness(nancycats(), "N100", "N111", alleles = allele_frequencies, inbreeding = true, verbose = true)

## Calculate theta and r
"""
    relatedness_dyadicML(Δ::Array{Float64,2})
Takes the Δ coefficents (with or without inbreeding allowed) and calculates the coefficient of relatedness
"""
function relatedness_from_Δ(Δ::Array{Float64,2})
    θ = Δ[1] + 0.5 * (Δ[3] + Δ[5] + Δ[7]) + 0.25 * Δ[8]
    2 * θ
end

#relatedness_dyadicML(Δ_inbreeding)
#relatedness_dyadicML(Δ_noinbreeding)
#Relatedness R package appears to have a bug. When allow.inbreeding = TRUE the relatedness value is the same as when I assume no inbreeding
#when you set allow.inbreeding = FALSE then the relatedness calculated is the same as when I assume there is inbreeding


"""
    dyadicML_relatedness(data::PopObj, ind1::String, ind2::String; alleles::Dict, inbreeding::Bool = true, verbose::Bool = true)
Calculates the dyadic maximum likelihood relatedness using all available loci following
Milligan 2002 dyadic relatedness - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1462494/pdf/12663552.pdf

Inbreeding can either be assumed not to occur (inbreeding = false) or be included in the calculation of the relatedness coefficient
Verbose controls the verbosity of the optimization process to find the Maximum likelihood Δ coefficents

"""
function dyadicML_relatedness(data::PopObj, ind1::String, ind2::String; alleles::Dict, inbreeding::Bool = true, verbose::Bool = true)

    Pr_L_S = all_loci_Pr_L_S(data, ind1, ind2, alleles)

    if !inbreeding
        no_inbreeding(Pr_L_S)
    end

    Δ,convergence = Δ_optim(Pr_L_S, verbose)

    return relatedness_from_Δ(Δ), Δ, convergence
end





"""
    pairwise_relatedness(data::PopObj, method::String, inbreeding::Bool = true, verbose::Bool = true)
Calculates various pairwise relatedness measures between all pairs of individuals based on the entire sample population
allele frequency

If the method is able to account for inbreeding in it's calculation then that option may be used

Currently only the Milligan 2002 Dyadic Maximum Likelihood relatedness estimator is implemented

"""
function pairwise_relatedness(data::PopObj; method::String, inbreeding::Bool = true, verbose::Bool = true)
    allele_frequencies = Dict()
    for locus in names(data.loci)
        allele_frequencies[String(locus)] = allele_freq(data.loci[:, locus])
    end

    #Add switch to slightly change output depending on relatednes metric (e.g. convergence doesn't make sense for Moments Estimators)
    output = DataFrame(ind1 = [], ind2 = [], relatedness = [], convergence = [])
    for (i,ind1) in enumerate(data.samples.name)
        i == length(data.samples.name) && break
        for ind2 in i+1:length(data.samples.name)

            #In here add switch for which relatedness metric to calculate
            dyad_out = dyadicML_relatedness(data, ind1, data.samples.name[ind2], alleles = allele_frequencies, inbreeding = inbreeding, verbose = verbose)
            append!(output,DataFrame(ind1 = ind1, ind2 = data.samples.name[ind2], relatedness = dyad_out[1], convergence = dyad_out[3]))
        end
        if verbose
            println("All pairs with ", ind1, " finished")
        end
    end

    return output
end




cat_rel = pairwise_relatedness(nancycats(), method = "dyadml", inbreeding = true, verbose = false)


test = filter(row -> row.convergence != :Optimal, cat_rel)
test2 = filter(row -> row.convergence != :Suboptimal, test)
test3 = filter(row -> row.convergence != :UserLimit, test2)

#maxit = 100 feastol = 1e-7 - 37.4% non-optimal; 49 total hit user limit
#maxit = 1000 feastol = 1e-7 - 26% non-optimal;  21 total hit user limit - more backtracking issues
#maxit = 100 feastol = 1e-6 - 30.7% non-optimal; 30 total hit user limit 76 infeasable
#maxit = 100 feastol = 1e-5 - 23.2% non-optimal; 0 total hit user limit  5728 infeasable - way fewer backtracking issues (only 2)
#maxit = 100 feastol = 5e-6 - 23.1% non-optimal; 0 total hit user limit 1570 infeasable - 7 failed backtracking
#maxit = 100 feastol = 1e-4 - % non-optimal;  total hit user limit  infeasable -  failed backtracking

## Introduce whacky fun if is infeasable to drop the feastol to 1e-8
#maxit = 100 feastol = 5e-6 - 21.9% non-optimal; 49 total hit user limit 0 infeasable - 7 failed backtracking

## SCSSolver defaults


6133/27966
76-106
