## Method 1

function get_uncondensed_state(ind1::T, ind2::T) where T <: Tuple
    i,j = ind1
    k,l = ind2

    δ = 1 * (i == j == k == l) +
    2 * ((i == j == k) & (i != l)) +
    3 * ((i == j == l) & (i != k)) +
    4 * ((i == k == l) & (i != j)) +
    5 * ((j == k == l) & (i != j)) +
    6 * ((i == j) & (k == l) & (i != k)) +
    7 * ((i == j) & (k != l) & (i != k) & (i != l)) +
    8 * ((k == l) & (i != j) & (i != k) & (j != k)) +
    9 * ((i == k) & (j == l) & (i != j) & (k != l)) +
    10 * ((i == k) & (j != l) & (i != j) & (k != l)) +
    11 * ((j == l) & (i != k) & (i != j) & (k != l)) +
    12 * ((i == l) & (j == k) & (i != j) & (k != l)) +
    13 * ((i == l) & (j != k) & (i != j) & (k != l)) +
    14 * ((i != l) & (j == k) & (i != j) & (k != l)) +
    15 * ((i != j) & (i != k) & (i != l) & (j != k) & (j != l) & (k != l))

    return δ
end


s1(p) = [p[1], p[1]^2, p[1]^2, p[1]^3, p[1]^2, p[1]^3, p[1]^2, p[1]^3, p[1]^4]

s2(p) = [0.0, 0.0, prod(p[[1,4]]), 2.0 * prod(p[[1,4]]) * p[1], 0.0, 0.0, 0.0, prod(p[[1,4]]) * p[1], 2.0 * prod(p[[1,4]]) * p[1]^2] #D3/S2

s3(p) = [0.0, 0.0, prod(p[[1,3]]), 2.0 * prod(p[[1,3]]) * p[1], 0.0, 0.0, 0.0, prod(p[[1,3]]) * p[1], 2.0 * prod(p[[1,3]]) * p[1]^2] #D3/S3

s4(p) = [0.0, 0.0, 0.0, 0.0, prod(p[[1,2]]), 2.0 * prod(p[[1,2]]) * p[1], 0.0, prod(p[[1,2]]) *p[1], 2.0 * prod(p[[1,2]]) * p[1]^2] #D5/S4

s5(p) = [0.0, 0.0, 0.0, 0.0, prod(p[[1,2]]), 2.0 * prod(p[[1,2]]) * p[2], 0.0, prod(p[[1,2]]) *p[2], 2.0 * prod(p[[1,2]]) * p[2]^2] #D5/S5

s6(p) = [0.0, prod(p[[1,3]]), 0.0, prod(p[[1,3]]) * p[3], 0.0, prod(p[[1,3]]) * p[1], 0.0, 0.0, prod(p[[1,3]])^2] #D2/S6

s7(p) = [0.0, 0.0, 0.0, 2.0 * prod(p[[1,3,4]]), 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,3,4]]) * p[1]] #D4/S7

s8(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2,3]]), 0.0, 0.0, 2.0 * prod(p[[1,2,3]]) * p[3]] #D6/S8

s9(p) =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2]]), prod(p[[1,2]]) * sum(p[[1,2]]), 4.0 * prod(p[[1,2]])^2] #D7/S9

s10(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,4]]), 4.0 * prod(p[[1,2,4]]) * p[1]] #D8/S10

s11(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,3]]), 4.0 * prod(p[[1,2,3]]) * p[2]] #D8/S11

s12(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2]]), prod(p[[1,2]]) * sum(p[[1,2]]), 4.0 * prod(p[[1,2]])^2] #D7/S12

s13(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,3]]), 4.0 * prod(p[[1,2,3]]) * p[1]] #D8/S13

s14(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,4]]), 4.0 * prod(p[[1,2,4]]) * p[2]] #D8/S14

s15(p) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * prod(p)] #D9/S15

δ = (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15)


function probability_states(ind1::T, ind2::T, alleles::Dict) where T <: Tuple
    i,j = ind1
    k,l = ind2
    p = [alleles[i], alleles[j], alleles[k], alleles[l]]

    return δ[get_uncondensed_state(ind1, ind2)](p)
end



s1(p) = (p[1], p[1]^2, p[1]^2, p[1]^3, p[1]^2, p[1]^3, p[1]^2, p[1]^3, p[1]^4)

s2(p) = (0.0, 0.0, prod(p[[1,4]]), 2.0 * prod(p[[1,4]]) * p[1], 0.0, 0.0, 0.0, prod(p[[1,4]]) * p[1], 2.0 * prod(p[[1,4]]) * p[1]^2) #D3/S2

s3(p) = (0.0, 0.0, prod(p[[1,3]]), 2.0 * prod(p[[1,3]]) * p[1], 0.0, 0.0, 0.0, prod(p[[1,3]]) * p[1], 2.0 * prod(p[[1,3]]) * p[1]^2) #D3/S3

s4(p) = (0.0, 0.0, 0.0, 0.0, prod(p[[1,2]]), 2.0 * prod(p[[1,2]]) * p[1], 0.0, prod(p[[1,2]]) *p[1], 2.0 * prod(p[[1,2]]) * p[1]^2) #D5/S4

s5(p) = (0.0, 0.0, 0.0, 0.0, prod(p[[1,2]]), 2.0 * prod(p[[1,2]]) * p[2], 0.0, prod(p[[1,2]]) *p[2], 2.0 * prod(p[[1,2]]) * p[2]^2) #D5/S5

s6(p) = (0.0, prod(p[[1,3]]), 0.0, prod(p[[1,3]]) * p[3], 0.0, prod(p[[1,3]]) * p[1], 0.0, 0.0, prod(p[[1,3]])^2) #D2/S6

s7(p) = (0.0, 0.0, 0.0, 2.0 * prod(p[[1,3,4]]), 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,3,4]]) * p[1]) #D4/S7

s8(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2,3]]), 0.0, 0.0, 2.0 * prod(p[[1,2,3]]) * p[3]) #D6/S8

s9(p) =  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2]]), prod(p[[1,2]]) * sum(p[[1,2]]), 4.0 * prod(p[[1,2]])^2) #D7/S9

s10(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,4]]), 4.0 * prod(p[[1,2,4]]) * p[1]) #D8/S10

s11(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,3]]), 4.0 * prod(p[[1,2,3]]) * p[2]) #D8/S11

s12(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p[[1,2]]), prod(p[[1,2]]) * sum(p[[1,2]]), 4.0 * prod(p[[1,2]])^2) #D7/S12

s13(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,3]]), 4.0 * prod(p[[1,2,3]]) * p[1]) #D8/S13

s14(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p[[1,2,4]]), 4.0 * prod(p[[1,2,4]]) * p[2]) #D8/S14

s15(p) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * prod(p)) #D9/S15

δ_t = (s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15)


function probability_states_tuple(ind1::T, ind2::T, alleles::Dict) where T <: Tuple
    i,j = ind1
    k,l = ind2
    p = [alleles[i], alleles[j], alleles[k], alleles[l]]

    return δ_t[get_uncondensed_state(ind1, ind2)](p)
end



## Method 2
S1 = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
S2 = [1 2; 3 4]
S3 = [1 2; 1 3; 2 3]
S4 = [1 2]
S5 = [1 3; 1 4; 3 4]
S6 = [3 4]
S7 = [1 3; 2 4]
S8 = [1 3]
S9 = [1 3][[1 3][:,1] .!= 1, :]
ibd_states = [S1, S2, S3, S4, S5, S6, S7, S8, S9]

function get_jacquard_state(ind1::T, ind2::T, ibd_states = ibd_states) where T <: Tuple
    the_shared = [intersect(ind1, ind2)...,0]
    the_lone = [setdiff(ind1, ind2)..., setdiff(ind2, ind1)...,0]

    ind1_resort = [ind1[(sum(ind1 .∈ the_shared',dims=2) .>= 1)[:,1]]..., ind1[(sum(ind1 .∈ the_lone',dims=2) .>= 1)[:,1]]...]
    ind2_resort = [ind2[(sum(ind2 .∈ the_shared',dims=2) .>= 1)[:,1]]..., ind2[(sum(ind2 .∈ the_lone',dims=2) .>= 1)[:,1]]...]

    adj_mat = Array{Int8}(undef, 6, 2)
    idx = 0
    for i in 1:3, j in i+1:4
        idx += 1
        adj_mat[idx,1:2] = [i,j] .* (1.0 * ([ind1_resort..., ind2_resort...][i] == [ind1_resort..., ind2_resort...][j]))
    end
    adj_mat = adj_mat[adj_mat[:,1] .> 0,:]

    jaq_state = 0
    for i in 1:9
        jaq_state += i * (adj_mat == ibd_states[i])
    end
    return jaq_state
end


## Single type?
struct JacquardPair
    genotype1::Genotype
    genotype2::Genotype
    state::Int
end

function probability_state(genos::JacquardPair, alleles::Dict)::Vector{Float64}
    state = genos.state
    if state == 1
        p = alleles[genos.genotype1[1]]
        [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]
    elseif state == 2
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]])
        [0.0, prod(p), 0.0, prod(p) * p[2], 0.0, prod(p) * p[1], 0.0, 0.0, prod(p)^2]

    elseif state == 3
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[2]])
        [0.0, 0.0, prod(p), 2.0 * prod(p) * p[1], 0.0, 0.0, 0.0, prod(p) * p[1], 2.0 * prod(p) * p[1]^2]

    elseif state == 4
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
        [0.0, 0.0, 0.0, 2.0 * prod(p), 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p) * p[1]]

    elseif state == 5
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
        [0.0, 0.0, 0.0, 0.0, prod(p), 2.0 * prod(p) * p[1], 0.0, prod(p) *p[1], 2.0 * prod(p) * p[1]^2]

    elseif state == 6
        p = (alleles[genos.genotype2[1]], alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
        [0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p), 0.0, 0.0, 2.0 * prod(p) * p[1]]

    elseif state == 7
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p), prod(p) * sum(p), 4.0 * prod(p)^2]

    elseif state == 8
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[2]])
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p), 4.0 * prod(p) * p[1]]

    elseif state == 9
        p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0 * prod(p)]
    else
        [-9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0]
    end
end


## multiple types?
struct JacquardZero
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardOne
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardTwo
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardThree
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardFour
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardFive
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardSix
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardSeven
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardEight
    genotype1::Genotype
    genotype2::Genotype
end

struct JacquardNine
    genotype1::Genotype
    genotype2::Genotype
end

function probability_state(genos::JacquardZero, alleles::Dict)::Vector{Float64}
    [-9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0, -9.0]
end

function probability_state(genos::JacquardOne, alleles::Dict)::Vector{Float64}
    p = alleles[genos.genotype1[1]]
    [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]
end

function probability_state(genos::JacquardTwo, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]])
    [0.0, prod(p), 0.0, prod(p) * p[2], 0.0, prod(p) * p[1], 0.0, 0.0, prod(p) * prod(p)]
end

function probability_state(genos::JacquardThree, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]])
    [0.0, 0.0, prod(p), 2.0 * prod(p) * p[1], 0.0, 0.0, 0.0, prod(p) * p[1], 2.0 * prod(p) * p[1]^2]
end

function probability_state(genos::JacquardFour, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
    [0.0, 0.0, 0.0, 2.0 * prod(p), 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p) * p[1]]
end

function probability_state(genos::JacquardFive, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[2]], alleles[genos.genotype1[1]])
    [0.0, 0.0, 0.0, 0.0, prod(p), 2.0 * prod(p) * p[1], 0.0, prod(p) *p[1], 2.0 * prod(p) * p[1]^2]
end

function probability_state(genos::JacquardSix, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype2[1]], alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
    [0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p), 0.0, 0.0, 2.0 * prod(p) * p[1]]
end

function probability_state(genos::JacquardSeven, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 * prod(p), prod(p) * sum(p), 4 * prod(p) * prod(p)]
end

function probability_state(genos::JacquardEight, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[2]])
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, prod(p), 4 * prod(p) * p[1]]
end

function probability_state(genos::JacquardNine, alleles::Dict)::Vector{Float64}
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4 * prod(p)]
end
