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
        @test missingdata(cats, by = "full") isa DataFrame
    end

    @testset "pairwise identical" begin
        pw_i = pairwiseidentical(cats)
        @test pw_i isa DataFrame
        @test size(pw_i) == (27966, 4)
        pw_i_2 = pairwiseidentical(cats, cats.sampleinfo.name[1:10])
        @test pw_i_2 isa DataFrame
        @test size(pw_i_2) == (45, 4)
    end

    @testset "frequency tables" begin
        @test genofreqtable(cats) isa DataFrame
        @test genofreqtable(cats, by = "population") isa DataFrame
        @test allelefreqtable(cats) isa DataFrame
        @test allelefreqtable(cats, by = "population") isa DataFrame
    end
end

end # module