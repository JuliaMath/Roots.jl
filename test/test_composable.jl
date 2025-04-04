# Simple tests to make sure Roots works with other
# packages. This is not part of runtests, as no expectation
# that CI should run these

using Roots
using Test
using Unitful
using Polynomials
using ForwardDiff
using Measurements

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
        g = (9 + 8 // 10) * m / s^2
        v0 = 10m / s
        y0 = 16m
        y(t) = -g * t^2 + v0 * t + y0

        for order in orders
            @test find_zero(y, 1.8s, order) ≈ 1.886053370668014s
            @test find_zero(y, 1.8f0s, order) isa typeof(1.88f0s)
        end

        for M in [Roots.Bisection(), Roots.A42(), Roots.AlefeldPotraShi()]
            @test find_zero(y, (1.8s, 1.9s), M) ≈ 1.886053370668014s
            @test find_zero(y, (1.8f0s, 1.9f0s), M) isa typeof(1.88f0s)
        end

        xrts = find_zeros(y, 0s, 10s)
        @test length(xrts) == 1
        @test xrts[1] ≈ 1.886053370668014s

        # issue #434
        xzs1 = find_zeros(x -> cos(x / 1u"m"), -1.6u"m", 2u"m")
        @test length(xzs1) == 2 && maximum(xzs1) ≈ 1.5707963267948966 * u"m"

        FX = ZeroProblem(y, (0.0f0s, 2.0f0s))
        prob = Roots.init(FX, Roots.AlefeldPotraShi())
        @test Roots.is_small_Δx(prob.M, prob.state, prob.options) isa Bool  # does not throw
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

    # Measurements # issue #453
    @testset "Measurements.jl" begin
        a = measurement(1.0, 0.1)
        b = measurement(-3.0, 0.1)
        c = measurement(-10.0, 0.1)
        f(x) = a * x^2 + b * x + c
        x₀ = (measurement(-3.0, 0.1), measurement(0.0, 0.1))
        for M in (A42(), AlefeldPotraShi(), Bisection())
            @test find_zero(f, x₀, M) ≈ -2.0
        end
        @test find_zero(f, measurement(0.0, 0.1)) ≈ -2.0
    end
end
