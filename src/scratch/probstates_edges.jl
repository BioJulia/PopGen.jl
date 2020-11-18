function probstates(ind1::Tuple, ind2::Tuple)
    ListOfTuples = [ind1, ind2]
    Edge = Dict{Tuple{Int64,Int64},Int64}()
    Ct = 1
    for i in 1:4
       for j in (i+1):4
          Edge[(i,j)] = Ct
          Ct += 1
       end
    end
    return Edge
    A = zeros(6)
    for (i, j) in ListOfTuples
        if i > j
            i, j = j, i
        end
        A[Edge[(i,j)]] = 1
    end
    A == [1,1,1,1,1,1] && return 1
    A == [1,0,0,0,0,1] && return 2
    A == [1,1,0,1,0,0] && return 3
    A == [1,0,0,0,0,0] && return 4
    A == [0,1,1,0,0,1] && return 5
    A == [0,0,0,0,0,1] && return 6
    A == [0,1,0,0,1,0] && return 7
    A == [0,1,0,0,0,0] && return 8
    A == [0,0,0,0,0,0] && return 9
end

if A[1] == 1
    if A[2] == 1
        A[5] == 1 ? return 1 : return 3
    else
        A[6] == 1 ? return 2 : return 4
    end
else
    if A[2] == 1
        A[3] == 1 && return 5
        A[5] == 1 ? return 7 : return 8
    else
        A[6] == 1 ? return 6 : return 9
    end
end


####################
####################
#=
Calculate Pr(Li | Sj)
If the allele identity falls into this class (L1-L9), generate the
probabilities of it belonging to each of the different classes and
return that array of 9 distinct probabilities
=#
function probability_state_table(x::Tuple, y::Tuple, alleles::Dict)
    #TODO Improve how groups are decided based on how similar things are done with moments estimators
    x1, x2 = x
    y1, y2 = y
    ## class L1 -  AᵢAᵢ AᵢAᵢ ##
    if x1 == x2 == y1 == y2
        p = alleles[x1]
        (p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4)

    ## class L2 - AᵢAᵢ AⱼAⱼ ##
    elseif x1 == x2 != y1 == y2
        # (x1 == x2) & (y1 == y2) & (x1 != y1)
        p1,p2 = alleles[x1], alleles[y1]
        prod_p = prod((p1,p2))
        [0, prod_p, 0, prod_p * p2, 0, prod_p * p1, 0, 0, prod_p^2]

    ## class L3a - AᵢAᵢ AᵢAⱼ ## - has issues because of allele order
    elseif x1 == x2 == y1 != y2
        #((x1 == x2 == y1) & (x1 != y2))
        p1,p2 = alleles[x1], alleles[y2]
        prod_p = prod(p)
        [0, 0, prod_p, 2 * prod_p * p1, 0, 0, 0, prod_p * p1, 2 * prod_p * p1^2]

    ## class L3b - AᵢAᵢ AⱼAᵢ ## - has issues because of allele order
    elseif x1 == x2 == y2 != y1
        p1,p2 = alleles[x1], alleles[y1]
        prod_p = prod(p)
        [0, 0, prod_p, 2 * prod_p * p1, 0, 0, 0, prod_p * p1, 2 * prod_p * p1^2]

    ## class L4 - AᵢAᵢ AⱼAₖ ##
    elseif x1 == x2 != y1 != y2 != x1
        #(x1 == x2) & (y1 != y2) & (x1 != y1) & (x1 != y2)
        p = (alleles[x1], alleles[y1], alleles[y2])
        prod_p = prod(p)
        [0, 0, 0, 2 * prod_p, 0, 0, 0, 0, 2 * prod_p * p[1]]

    ## L5a - AiAj AiAi ## - has issues because of allele order
    elseif x1 == y1 == y2 != x2
        p1,p2 = alleles[x1], alleles[x2]
        prod_p = prod(p)
        [0, 0, 0, 0, prod_p, 2 * prod_p * p1, 0, prod_p * p1, 2 * prod_p * p1^2]

    ## L5b - AjAi AiAi ## - has issues because of allele order
    elseif x2 == y1 == y2 != x2
        p1,p2 = alleles[x2], alleles[x1]
        prod_p = prod(p)
        [0, 0, 0, 0, prod_p, 2 * prod_p * p1, 0, prod_p *p1, 2 * prod_p * p1^2]

    ## L6 - AjAk AiAi ##
    elseif y1 != x2 != x1 != y1 == y2
        #(x1 != x2) & (y1 == y2) & (x1 != y1) & (x2 != y1)
        p = (alleles[y1], alleles[x1], alleles[x2])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 2 * prod_p, 0, 0, 2 * prod_p * p[1]]

    ## L7 - AiAj AiAj ##
    elseif y1 == x1 != x2 == y2 
        #(x1 == y1) & (x2 == y2) & (x1 != x2)
        p = (alleles[x1], alleles[x2])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 0, 2 * prod_p, prod_p * sum(p), 4 * prod_p * prod_p]

    ## L8a - AiAj AiAk ##  - has issues because of allele order
    elseif x2 != x1 == y1 != y2 != x2
        #(x1 == y1) & (x1 != x2) & (y1 != y2) & (x2 != y2)
        p = (alleles[x1], alleles[x2], alleles[y2])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 0, 0, prod_p, 4 * prod_p * p[1]]

    ## L8b - AjAi AkAi ##  - has issues because of allele order
    elseif x1 != x2 == y2 != y1 != x1
        #(x2 == y2) & (x1 != x2) & (y1 != y2) & (x1 != y1)
        p = (alleles[x2], alleles[x1], alleles[y1])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 0, 0, prod_p, 4 * prod_p * p[1]]

    ## L8c - AjAi AiAk ##  - has issues because of allele order
    elseif x1 != x2 == y1 != y2 != x1 
        #(x2 == y1) & (x1 != x2) & (y1 != y2) & (x1 != y2)
        p = (alleles[x2], alleles[x1], alleles[y2])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 0, 0, prod_p, 4 * prod_p * p[1]]

    ## L8d - AiAj AkAi ##  - has issues because of allele order
    elseif x2 != x1 == y2 != y2 & (x1 != y1)
        #(x1 == y2) & (x1 != x2) & (y1 != y2) & (x1 != y1)
        p = (alleles[x1], alleles[x2], alleles[y1])
        prod_p = prod(p)
        [0, 0, 0, 0, 0, 0, 0, prod_p, 4 * prod_p * p[1]]

    ## L9 - AiAj AkAl ##
    elseif x1 != x2 != y1 != y2 != x1 != y1 & (y2 != x2)
        #(x1 != x2) & (x1 != y1) & (x1 != y2) & (x2 != y1) & (x2 != y2) & (y1 != x2)
        p = (alleles[x1], alleles[x2], alleles[y1], alleles[y2])
        [0, 0, 0, 0, 0, 0, 0, 0, 4 * prod(p)]
    else
        [-9, -9, -9, -9, -9, -9, -9, -9, -9]
    end
end