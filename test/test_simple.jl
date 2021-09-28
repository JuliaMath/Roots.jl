using Test
using BenchmarkTools

@testset "simpler implementations" begin

    # bisection
    xrt = Roots.bisection(sin, 3.0, 4.0)
    @test isapprox(xrt, pi)

    xrt = Roots.bisection(sin, 3.0, 4.0, xatol=1e-3)
    @test abs(sin(xrt)) >= 1e-7  # not to0 close

    xrt = Roots.bisection(sin, big(3.0), big(4.0))
    @test isapprox(xrt, pi)

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

    # dfree
    fpoly = x -> x^5 - x - 1
    xrt = Roots.dfree(fpoly, 1.0)
    @test abs(fpoly(xrt)) <= 1e-14

    # newton
    @test Roots.newton((sin, cos), 3.0) ≈ pi
    u = Roots.newton(x -> (sin(x), sin(x) / cos(x)), 3.0, xatol=1e-10, xrtol=1e-10)
    @test abs(u - pi) <= 1e-8
end
