S1 = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]
S2 = [1 2; 3 4]
S3 = [1 2; 1 3; 2 3]
S4 = [1 2]
S5 = [1 3; 1 4; 3 4]
S6 = [3 4]
S7 = [1 3; 2 4]
S8 = [1 3]
S9 = [1 3][[1 3][:,1] .!= 1, :]

function get_jacquard_state(ind1::T, ind2::T) where T <: Tuple
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
    return adj_mat[adj_mat[:,1] .> 0,:]
end


## Single type?
struct JacquardPair
    genotype1::Genotype
    genotype2::Genotype
    state::Int
end

function probability_state(genos::JacquardPair, state::genos.state, alleles::Dict)::Vector{Float64}
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
