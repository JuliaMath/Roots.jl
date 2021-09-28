# Simple tests to make sure Roots works with other
# packages. This is not part of runtests, as no expectation
# that CI should run these

using Roots
using Test
using Unitful
using Polynomials
using ForwardDiff

@testset "Test composability with other packages" begin
    orders = [
        Order0(),
        Order1(),
        Roots.Order1B(),
        Order2(),
        Roots.Order2B(),
        Order5(),
        Order8(),
        Order16(),
    ]

    # Unitful
    @testset "find zero(s) with Unitful" begin
        s = u"s"
        m = u"m"
        g = 9.8 * m / s^2
        v0 = 10m / s
        y0 = 16m
        y(t) = -g * t^2 + v0 * t + y0

        for order in orders
            @test find_zero(y, 1.8s, order) ≈ 1.886053370668014s
        end

        for M in [Roots.Bisection(), Roots.A42(), Roots.AlefeldPotraShi()]
            @test find_zero(y, (1.8s, 1.9s), M) ≈ 1.886053370668014s
        end

        xrts = find_zeros(y, 0s, 10s)
        @test length(xrts) == 1
        @test xrts[1] ≈ 1.886053370668014s
    end

    # Polynomials
    @testset "find zero(s) with Polynomials" begin
        m = s = 1.0
        g = 9.8 * m / s^2
        v0 = 10m / s
        y0 = 16m
        y(t) = -g * t^2 + v0 * t + y0

        x = variable()
        for order in orders
            @test find_zero(y(x), 1.8, order) ≈ 1.8860533706680143
        end

        for M in [Roots.Bisection(), Roots.A42(), Roots.AlefeldPotraShi()]
            @test find_zero(y(x), (1.8, 1.9), M) ≈ 1.8860533706680143
        end

        xrts = find_zeros(y(x), 0, 10)
        @test length(xrts) == 1
        @test xrts[1] ≈ 1.8860533706680143
    end

    # ForwardDiff
    # taking a derivative of a function that using find_zero.
    D(f) = x -> ForwardDiff.derivative(f, x)
    F(x, y) = (x + 1) * (x^2 + y^2) - 2x^2 # has  loop for x ∈ (0,1)
    h(x) = find_zero(y -> F(x, y), -1 / 4 * one(x))  # one(x) to get proper type for D
    α = find_zero(D(h), (0, 1)) # find lowest point on loop
    @test h(α) ≤ h(α + 0.1)
    @test h(α) ≤ h(α - 0.1)
end
