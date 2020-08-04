
using Convex
using ECOS


"""
    probability_state_table(x::Tuple, y::Tuple, loc::Symbol, alleles::Dict)
Calculate the probability of observing the particular allele state given each of
the 9 Jacquard Identity States for a single locus to create Table 1 from
Milligan 2003.
"""
function probability_state_table(x::Tuple, y::Tuple, alleles::Dict)
    #TODO Improve how groups are decided based on how similar things are done with moments estimators
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
    Δ_optim(Pr_L_S::Transpose{Float64,Array{Float64,2}}, verbose::Bool)
Takes the probability of the allelic state given the identity by descent from
all available loci (either allowing for inbreeding or not) and calculated the
maximum likelihood Δ coefficients
"""
function Δ_optim(Pr_L_S::Array{Float64,2}, verbose::Bool = false)
    #Δ is what needs to be optimized
    #consist of 9 values between 0 and 1 which must also sum to 1
    #is then used to calculate relatedness

    Δ = Variable(9)
    problem = maximize(sum(log(Pr_L_S * Δ)))
    problem.constraints += 0 <= Δ[1:9]
    problem.constraints += Δ[1:9] <= 1
    problem.constraints += sum(Δ) <= 1
    problem.constraints += 0 <= sum(Δ)
    problem.constraints += 2 * (Δ[1] + 0.5 * (Δ[3] + Δ[5] + Δ[7]) + 0.25 * Δ[8]) <= 1
    problem.constraints += 0 <= 2 * (Δ[1] + 0.5 * (Δ[3] + Δ[5] + Δ[7]) + 0.25 * Δ[8])

    Convex.solve!(problem, ECOS.Optimizer(verbose = verbose, maxit=100), verbose = verbose) #maxit=100,
    #Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit=100, feastol=5e-6, reltol = 1e-3, reltol_inacc = 5e-2), verbose = verbose)
    #Convex.solve!(problem, SCSSolver(verbose = verbose, max_iters = 100), verbose = verbose)

    Δ.value, problem.status
    # Should probably include some output that confirms that it did in fact
    # converge and/or use multiple random starts to confirm not a local maxima
end

#### No inbreeding assumption
function Δ_optim_noInbreeding(Pr_L_S::Array{Float64,2}, verbose::Bool = false)
    #Δ is what needs to be optimized
    #consist of 9 values between 0 and 1 which must also sum to 1
    #is then used to calculate relatedness

    Δ = Variable(3)
    problem = maximize(sum(log(Pr_L_S * Δ)))
    problem.constraints += 0 <= Δ[1:3]
    problem.constraints += Δ[1:3] <= 1
    problem.constraints += sum(Δ) <= 1
    problem.constraints += 0 <= sum(Δ)
    problem.constraints += 2 * (0.5 * Δ[1] + 0.25 * Δ[2]) <= 1
    problem.constraints += 0 <= 2 * (0.5 * Δ[1] + 0.25 * Δ[2])

    Convex.solve!(problem, ECOS.Optimizer(verbose = verbose, maxit=100), verbose = verbose) #maxit=100,
    #Convex.solve!(problem, ECOSSolver(verbose = verbose, maxit=100, feastol=5e-6, reltol = 1e-3, reltol_inacc = 5e-2), verbose = verbose)
    #Convex.solve!(problem, SCSSolver(verbose = verbose, max_iters = 100), verbose = verbose)

    Δ.value, problem.status
    # Should probably include some output that confirms that it did in fact
    # converge and/or use multiple random starts to confirm not a local maxima
end


"""
dyadicLikelihood(ind1::GenoArray, ind2::GenoArray, locus_names::Vector{Symbol}; alleles::NamedTuple)
Calculates the maximum likelihood based relatedness using all available loci following following Milligan (2002)

-Single Locus Equation:
-How to combine multiple loci: NA inherently multi-locus

-Assumes inbreeding can be present

Milligan, B. G. (2003). Maximum-likelihood estimation of relatedness. Genetics, 163(3), 1153-1167.
"""
function dyadicLikelihood(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    kw_dict = Dict(kwargs...)
    Pr_Ls = Array{Float64}(undef, length(locus_names), 9)

    idx = 0
    @inbounds for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        idx += 1
        Pr_Ls[idx, :] = probability_state_table(gen1, gen2, alleles[loc])
    end

    if kw_dict[:inbreeding] == true
        Δ = Δ_optim(Pr_Ls)
        θ = Δ[1][1] + 0.5 * (Δ[1][3] + Δ[1][5] + Δ[1][7]) + 0.25 * Δ[1][8]
    else
        Pr_Ls = Pr_Ls[:, 7:9]
        Δ = Δ_optim_noInbreeding(Pr_Ls)
        θ = (0.5 * Δ[1][1] + 0.25 * Δ[1][2])
    end
    return 2 * θ#, Δ[1], Δ[2] #secondary outputs for convergence & Δvalues if we want
end



#= should be redundant now
function dyadicLikelihood_noInbreeding(ind1::T, ind2::T, locus_names::Vector{Symbol}, alleles::U; kwargs...) where T <: GenoArray where U <: NamedTuple
    #TODO Add inbreeding toggle
    Pr_Ls = Array{Float64}(undef, length(locus_names), 9)

    idx = 0
    for (loc,gen1,gen2) in zip(locus_names, ind1, ind2)
        idx += 1
        Pr_Ls[idx, :] = probability_state_table(gen1, gen2, allele_frequencies[loc])
    end
    Pr_Ls = Pr_Ls[:, 7:9]

    Δ = Δ_optim_noInbreeding(Pr_Ls)

    θ = (0.5 * Δ[1][1] + 0.25 * Δ[1][2])
    return 2 * θ#, Δ[1], Δ[2] #secondary outputs for convergence & Δvalues if we want
end
=#