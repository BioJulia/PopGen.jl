
function JacquardOne
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

function probability_state(genos::JacquardOne, alleles::Dict)
    p = alleles[genos.genotype1[1]]
    [p, p^2, p^2, p^3, p^2, p^3, p^2, p^3, p^4]
end

function probability_state(genos::JacquardTwo, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]])
    [0, prod(p), 0, prod(p) * p[2], 0, prod(p) * p[1], 0, 0, prod(p) * prod(p)]
end

function probability_state(genos::JacquardThree, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]])
    [0, 0, prod(p), 2 * prod(p) * p[1], 0, 0, 0, prod(p) * p[1], 2 * prod(p) * p[1]^2]
end

function probability_state(genos::JacquardFour, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
    [0, 0, 0, 2 * prod(p), 0, 0, 0, 0, 2 * prod(p) * p[1]]
end

function probability_state(genos::JacquardFive, alleles::Dict)
    p = (alleles[genos.genotype1[2]], alleles[genos.genotype1[1]])
    [0, 0, 0, 0, prod(p), 2 * prod(p) * p[1], 0, prod(p) *p[1], 2 * prod(p) * p[1]^2]
end

function probability_state(genos::JacquardSix, alleles::Dict)
    p = (alleles[genos.genotype2[1]], alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
    [0, 0, 0, 0, 0, 2 * prod(p), 0, 0, 2 * prod(p) * p[1]]
end

function probability_state(genos::JacquardSeven, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]])
    [0, 0, 0, 0, 0, 0, 2 * prod(p), prod(p) * sum(p), 4 * prod(p) * prod(p)]
end

function probability_state(genos::JacquardEight, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[2]])
    [0, 0, 0, 0, 0, 0, 0, prod(p), 4 * prod(p) * p[1]]
end

function probability_state(genos::JacquardNine, alleles::Dict)
    p = (alleles[genos.genotype1[1]], alleles[genos.genotype1[2]], alleles[genos.genotype2[1]], alleles[genos.genotype2[2]])
    [0, 0, 0, 0, 0, 0, 0, 0, 4 * prod(p)]
end