using Test
using BenchmarkTools
import Roots: newton, halley, superhalley, quadratic_inverse, chebyshev_like

@testset "simpler implementations" begin

    # secant_method
    fpoly = x -> x^5 - x - 1
    xrt = Roots.secant_method(fpoly, 1.0)
    @test abs(fpoly(xrt)) <= 1e-15

    xrt = Roots.secant_method(fpoly, (1, 2))
    @test abs(fpoly(xrt)) <= 1e-14

    # muller
    fpoly = x -> x^5 - x - 1
    xrt = Roots.muller(fpoly, 1.0)
    @test xrt isa Real
    @test abs(fpoly(xrt)) <= 1e-15

    @test_throws DomainError Roots.muller(fpoly, -1.0)
    xrt = Roots.muller(fpoly, -1.0 + 0im)
    @test xrt isa Complex
    @test abs(fpoly(xrt)) <= 1e-15

    @test Roots.muller(cos, 1.0) ≈ π / 2
    expoly(z) = log(-z) * asin(z) / tanh(z)

    @test Roots.muller(expoly, -0.7 - 0.5im) ≈ -1.0

    # newton/halley/superhalley/quadratic_inverse/chebyshev_like
    @test Roots.newton((sin, cos), 3.0) ≈ pi
    u = Roots.newton(x -> (sin(x), sin(x) / cos(x)), 3.0, xatol=1e-10, xrtol=1e-10)
    @test abs(u - pi) <= 1e-8

    @test abs(newton(sin, cos, 0.5) - 0.0) <= 100 * eps(1.0)
    @test newton(cos, x -> -sin(x), 1.0) ≈ pi / 2
    @test newton(x -> x^2 - 2x - 1, x -> 2x - 2, 3.0) ≈ 2.414213562373095
    @test abs(newton(x -> exp(x) - cos(x), x -> exp(x) + sin(x), 3.0) - 0.0) <= 1e-14
    @test halley(x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2, 3.0) ≈ 2.414213562373095
    @test quadratic_inverse(x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2, 3.0) ≈ 2.414213562373095
    @test superhalley(x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2, 3.0) ≈ 2.414213562373095
    @test chebyshev_like(x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2, 3.0) ≈ 2.414213562373095
    a = halley(x -> exp(x) - cos(x), x -> exp(x) + sin(x), x -> exp(x) + cos(x), 3.0)
    @test abs(a - 0.0) <= 1e-14
    @test_throws Roots.ConvergenceFailed Roots.newton((x -> x^2 + 1, x -> 2x), 0)

    ## test with Complex input
    @test real(Roots.newton(x -> x^3 - 1, x -> 3x^2, 1 + im)) ≈ 1.0
    @test real(Roots.newton(x -> x^3 - 1, x -> 3x^2, 1 + 10im)) ≈ (-1 / 2)

    ## Issue #143 test with new interface
    Roots.newton(sin, cos, 3.0) ≈ π # uses find_zero
    Roots.newton((sin, cos), 3.0) ≈ π # uses simple

    fdf = x -> (sin(x), sin(x) / cos(x))  # (f, f/f')
    @test Roots.find_zero(fdf, 3.0, Roots.Newton()) ≈ π # uses find_zero
    Roots.newton(fdf, 3.0) ≈ π # uses simple




    # a42
    # simple a42()
    #m = run_tests(Roots.a42) # in test_bracketing
    #VERSION >= v"1.6" && @test isempty(m.failures)
    #@test m.evalcount <= 3000 # paper says 2884, this has 2877

end


@testset "simple: zero allocations" begin
    @test BenchmarkTools.@ballocated(Roots.secant_method(sin, 3)) == 0
    @test BenchmarkTools.@ballocated(Roots.muller(sin, 2.9, 3.0, 3.1)) == 0
    @test BenchmarkTools.@ballocated(Roots.newton((sin, cos), 3)) == 0
end
