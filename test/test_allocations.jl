using Test
import BenchmarkTools

@testset "solve: zero allocations" begin

    fs = (sin, cos, x -> -sin(x))
    x0 = (3,4)
    Ms = (Order1(), Order2(), Order5(), Order8(), Order16(),
          Roots.Order1B(), Roots.Order2B(),
          Roots.Bisection(), Roots.A42(), Roots.AlefeldPotraShi(), Roots.Brent(),
          Roots.Newton(), Roots.Halley(), Roots.Schroder()) # not Order0(), FalsePosition()

    for M ∈ Ms
        @test BenchmarkTools.@ballocated(solve(ZeroProblem($fs, $x0), $M)) == 0
    end


end

@testset "simple: zero allocations" begin
    @test BenchmarkTools.@ballocated( Roots.bisection(sin, 3,4) ) == 0
    @test BenchmarkTools.@ballocated( Roots.secant_method(sin, 3) ) == 0
    @test BenchmarkTools.@ballocated( Roots.muller(sin, 2.9,3.0,3.1) ) == 0
    @test BenchmarkTools.@ballocated( Roots.newton((sin,cos), 3) ) == 0
    @test BenchmarkTools.@ballocated( Roots.dfree(sin, 3) ) == 0
end
