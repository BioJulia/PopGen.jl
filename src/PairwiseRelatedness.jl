#= Steps to finish:

    Solve dyadic optimization issues
        combined backtracking failure - need to isolate and find cause
            #N111 N83
            #N63 N64
        ~1/3 - 1/5 not with optimal solution in cats dataset
    Multithread

=#

#= Would be good to include

    Implement alternative relatedness metrics

    Warning if a not implemented (or typo) of method included

    Streamline output

=#

"""
    pr_l_s(x::Tuple, y::Tuple, alleles::Dict)
Calculate the probability of observing the particular allele state given each of
the 9 Jacquard Identity States for a single locus to create Table 1 from
Milligan 2003.
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
    elseif (x[2] == y[1] == y[2] & (x[1] != x[2]))
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

    Δ = Variable(8)
    problem = maximize(sum(log(Pr_L_S * vcat(1 - sum(Δ), Δ))))
    problem.constraints += 0 <= 1 - sum(Δ)
    problem.constraints += 1 - sum(Δ) <= 1
    problem.constraints += 0 <= Δ[1:8]
    problem.constraints += Δ[1:8] <= 1

    #shifted from actual relatedness calculations because the 1 - sum(Δ) goes at beginning
    problem.constraints += 2 * ((1 - sum(Δ)) + 0.5 * (Δ[2] + Δ[4] + Δ[6]) + 0.25 * Δ[7]) <= 1
    problem.constraints += 0 <= 2 * ((1 - sum(Δ)) + 0.5 * (Δ[2] + Δ[4] + Δ[6]) + 0.25 * Δ[7])

    Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit = 100), verbose = verbose) #maxit=100,
    #Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit=100, feastol=5e-6, reltol = 1e-3, reltol_inacc = 5e-2), verbose = verbose)
    #Convex.solve!(problem, SCSSolver(verbose = verbose, max_iters = 100), verbose = verbose)

    vcat(1-sum(Δ.value), Δ.value), problem.status
    # Should probably include some output that confirms that it did in fact
    # converge and/or use multiple random starts to confirm not a local maxima
end

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

Inbreeding can either be assumed not to occur (inbreeding = false) or be
included in the calculation of the relatedness coefficient Verbose controls the
verbosity of the optimization process to find the Maximum likelihood Δ coefficents
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
    qg_relatedness(data::PopObj, ind1::String, ind2::String; alleles::Dict)
Calculates the moments based estimator of pairwise relatedness developed by Queller & Goodnight (1989).

Bases allele frequencies on entire population

Inbreeding can only be assumed not to exist.

See equation 3 in: https://www.nature.com/articles/hdy201752 for variant of estimator used
"""
function qg_relatedness(data::PopObj, ind1::String, ind2::String; alleles::Dict)

    n1 = n2 = d1 = d2 = 0

    for locus in String.(names(data.loci))
        #Extract the pair of interest's genotypes
        gen1 = get_genotype(data, sample = ind1, locus = locus)
        gen2 = get_genotype(data, sample = ind2, locus = locus)

        #Skip missing
        if gen1 !== missing && gen2 !== missing
            a = gen1[1]
            b = gen1[2]
            c = gen2[1]
            d = gen2[2]

            n1 = n1 + (Int(a==c) + Int(a==d) + Int(b==c) + Int(b==d) - 2 * (alleles[locus][a] + alleles[locus][b]))
            n2 = n2 + (Int(a==c) + Int(a==d) + Int(b==c) + Int(b==d) - 2 * (alleles[locus][c] + alleles[locus][d]))

            d1 = d1 + (2 * (1 + Int(a==b) - alleles[locus][a] - alleles[locus][b]))
            d2 = d2 + (2 * (1 + Int(c==d) - alleles[locus][c] - alleles[locus][d]))
        end
    end
    return (n1/d1 + n2/d2)/2

end



"""
    pairwise_relatedness(data::PopObj, method::String, inbreeding::Bool = true, verbose::Bool = true)
Calculates various pairwise relatedness measures between all pairs of individuals based on the entire sample population
allele frequency

If verbose is set to false then there is a progress bar. If set to true then there is estimator specific feedback and statements when an individual has been compared to all other pairs

If the method is able to account for inbreeding in it's calculation then that option may be used

Currently implemented are Milligan 2002 Dyadic Maximum Likelihood relatedness estimator and Queller & Goodnight 1989

"""
function pairwise_relatedness(data::PopObj; method::String, inbreeding::Bool = true, verbose::Bool = true)
    allele_frequencies = Dict()
    for locus in names(data.loci)
        allele_frequencies[String(locus)] = allele_freq(data.loci[:, locus])
    end

    if !verbose
        n = size(data.samples)[1]
        n = n*(n-1) ÷ 2
        prog = Progress(n, 1)
    end

    #Add switch to slightly change output depending on relatednes metric (e.g. convergence doesn't make sense for Moments Estimators)
    if method == "dyadml"
        output = DataFrame(ind1 = [], ind2 = [], relatedness = [], convergence = [])
    end
    if method == "qg"
        output = DataFrame(ind1 = [], ind2 = [], relatedness = [])
    end

    for (i,ind1) in enumerate(data.samples.name)
        i == length(data.samples.name) && break
        for ind2 in i+1:length(data.samples.name)

            if method == "dyadml"
                dyad_out = dyadicML_relatedness(data, ind1, data.samples.name[ind2], alleles = allele_frequencies, inbreeding = inbreeding, verbose = verbose)
                append!(output,DataFrame(ind1 = ind1, ind2 = data.samples.name[ind2], relatedness = dyad_out[1], convergence = dyad_out[3]))
            end
            if method == "qg"
                qg_out = qg_relatedness(data, ind1, data.samples.name[ind2], alleles = allele_frequencies)
                append!(output,DataFrame(ind1 = ind1, ind2 = data.samples.name[ind2], relatedness = qg_out))
            end

            if !verbose
                next!(prog)
            end
        end
        if verbose
            println("All pairs with ", ind1, " finished")
        end
    end

    return output
end

#=
cat_rel_noInbreeding = pairwise_relatedness(nancycats(), method = "dyadml", inbreeding = false, verbose = false)
cat_rel_Inbreeding = pairwise_relatedness(nancycats(), method = "dyadml", inbreeding = true, verbose = false)


cat_rel_qg = pairwise_relatedness(nancycats(), method = "qg", verbose = false)
test
=#
