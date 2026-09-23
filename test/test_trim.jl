using Roots
using Test

@testset "trim smoke test" begin
    f(x) = x^2 - 2
    @test find_zero(f, (1.0, 2.0)) ≈ sqrt(2)
    @test find_zero(sin, (3.0, 4.0)) ≈ pi
    @test solve(ZeroProblem(f, (1.0, 2.0)), Roots.Bisection()) ≈ sqrt(2)
end
