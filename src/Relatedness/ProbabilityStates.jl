## Single type?
struct JaqcuardPair
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