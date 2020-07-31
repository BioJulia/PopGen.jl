function get_jaquard_state(ind1::T, ind2::T) where T <: Tuple
    the_shared = [intersect(ind1, ind2)...,0]
    the_lone = [setdiff(ind1, ind2)..., setdiff(ind2, ind1)...,0]

    ind1_resort = [ind1[(sum(ind1 .∈ the_shared',dims=2) .>= 1)[:,1]]..., ind1[(sum(ind1 .∈ the_lone',dims=2) .>= 1)[:,1]]...]
    ind2_resort = [ind2[(sum(ind2 .∈ the_shared',dims=2) .>= 1)[:,1]]..., ind2[(sum(ind2 .∈ the_lone',dims=2) .>= 1)[:,1]]...]

    adj_mat = Array{Int8}(undef, 6, 2)
    idx = 0
    for i in collect(1:3)
        for j in (i + 1):4
            idx += 1
            adj_mat[idx,1:2] = [i,j] .* (1.0 * ([ind1_resort..., ind2_resort...][i] == [ind1_resort..., ind2_resort...][j]))
        end
    end
    return adj_mat[adj_mat[:,1] .> 0,:]
end




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
