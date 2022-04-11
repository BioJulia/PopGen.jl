module  TestAlleleMatrices

using PopGen
using Test

cats = @nancycats ;

@testset "AlleleMatrices.jl" begin
    @testset "_allelematrix" begin
        @test size(PopGen._allelematrix(cats, by = "count", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "missing", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "zero", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, by = "frequency", missings = "mean", scale = false, center = false)) == (237, 108)
        @test size(PopGen._allelematrix(cats, scale = true)) == (237, 108)
        @test size(PopGen._allelematrix(cats, center = true)) == (237, 108)
        @test PopGen._allelematrix(cats, scale = true) != PopGen._allelematrix(cats, center = true)
        a = PopGen._allelematrix(cats, by = "count")
        b = PopGen._allelematrix(cats, missings = "missing")
        c = PopGen._allelematrix(cats, missings = "zero")
        d = PopGen._allelematrix(cats)
        @test all(a .== b) == false
        @test all(a .== c) == false
        @test all(a .== d) == false
        @test all(b .== c) == false
        @test all(b .== d) == false
        @test all(c .== d) == false
    end
end

end # module