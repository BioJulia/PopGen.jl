module  TestDataExploration

using PopGen
using DataFrames
using Test

cats = @nancycats;

@testset "Data Exploration" begin

    @testset "missing data" begin
        @test missingdata(cats) isa DataFrame
        @test missingdata(cats, by = "population") isa DataFrame
        @test missingdata(cats, by = "locus") isa DataFrame
        @test missingdata(cats, by = "locusxpopulation") isa DataFrame
    end

    @testset "pairwise identical" begin
        pw_i = pairwiseidentical(cats)
        @test pw_i isa AbstractMatrix
        @test size(pw_i) == (237, 237)
        pw_i_2 = pairwiseidentical(cats, cats.sampleinfo.name[1:10])
        @test pw_i_2 isa AbstractMatrix
        @test size(pw_i_2) == (10, 10)
    end

    @testset "frequency tables" begin
        @test genofreqtable(cats) isa DataFrame
        @test genofreqtable(cats, by = "population") isa DataFrame
        @test allelefreqtable(cats) isa DataFrame
        @test allelefreqtable(cats, by = "population") isa DataFrame
    end
end

end # module