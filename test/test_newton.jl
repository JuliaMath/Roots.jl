using Test
import Roots.newton,
    Roots.halley, Roots.superhalley, Roots.quadratic_inverse, Roots.chebyshev_like

@testset "Test Newton, Halley, Schroder methods" begin
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

    # find_zero calls
    @test find_zero((x -> x^2 - 2x - 1, x -> 2x - 2), 3.0, Roots.Newton()) ≈
          2.414213562373095
    @test find_zero((x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2), 3.0, Roots.Halley()) ≈
          2.414213562373095
    @test find_zero((x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2), 3.0, Roots.SuperHalley()) ≈
          2.414213562373095
    @test find_zero((x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2), 3.0, Roots.ChebyshevLike()) ≈
          2.414213562373095
    @test find_zero(
        (x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2),
        3.0,
        Roots.QuadraticInverse(),
    ) ≈ 2.414213562373095
    @test find_zero((x -> x^2 - 2x - 1, x -> 2x - 2, x -> 2), 3.0, Roots.Schroder()) ≈
          2.414213562373095

    ## test with Complex input

    @test real(Roots.newton(x -> x^3 - 1, x -> 3x^2, 1 + im)) ≈ 1.0
    @test real(Roots.newton(x -> x^3 - 1, x -> 3x^2, 1 + 10im)) ≈ (-1 / 2)

    ## Issue #143 test with new interface
    Roots.newton(sin, cos, 3.0) ≈ π # uses find_zero
    Roots.newton((sin, cos), 3.0) ≈ π # uses simple

    fdf = x -> (sin(x), sin(x) / cos(x))  # (f, f/f')
    @test Roots.find_zero(fdf, 3.0, Roots.Newton()) ≈ π # uses find_zero
    Roots.newton(fdf, 3.0) ≈ π # uses simple

    fdfdf = x -> (sin(x), sin(x) / cos(x), -cos(x) / sin(x), sin(x) / cos(x))   # (f, f/f', f'/f'', f''/f''')
    @test Roots.find_zero(fdfdf, 3.0, Roots.Halley()) ≈ π

    # check that functions with multiple return values can work with other
    # methods
    for M in [
        Roots.Schroder(),
        Roots.Halley(),
        Roots.Newton(),
        Roots.Order1(),
        Roots.Order0(),
        Roots.QuadraticInverse(),
        Roots.SuperHalley(),
        Roots.ChebyshevLike(),
    ]
        @test Roots.find_zero(fdfdf, 3.0, M) ≈ π  # can pass function to others
    end

    for M in [Roots.Bisection(), Roots.A42(), Roots.AlefeldPotraShi()]
        @test Roots.find_zero(fdfdf, (3.0, 4.0), M) ≈ π  # can pass function to others
    end

    @test find_zero(
        x -> (x^2 - 2, (x^2 - 2) / 2x),
        1.0,
        Roots.Newton(),
        Roots.Bisection(),
    ) ≈ sqrt(2)

    @test_throws Roots.ConvergenceFailed Roots.newton((x -> x^2 + 1, x -> 2x), 0)
end

@testset "bracketed Halley" begin
    x₀, x̃₀, α = 1.0, 1.1, 1.1673039782614187
    f(x) = x^5 - x - 1
    fp(x) = 5x^4 - 1
    fpp(x) = 20x^3

    for M in (Roots.BracketedHalley(), Roots.BracketedChebyshev())
        @test find_zero((f, fp, fpp), (1, 2), M) ≈ α
    end
end

@testset "Lith Boonkkamp IJzerman methods" begin
    x₀, x̃₀, α = 1.0, 1.1, 1.1673039782614187
    f(x, p=1)  = x^5 - x - p
    fp(x)      = 5x^4 - 1
    fpp(x)     = 20x^3
    fppp(x)    = 60x
    fpppp(x)   = 60

    function fdf1(x, p=1)
        fx = f(x, p)
        ∂fx = fp(x)
        return fx, fx / ∂fx
    end

    function fdf2(x, p=1)
        fx = f(x, p)
        ∂fx = fp(x)
        ∂2fx = fpp(x)
        return fx, fx / ∂fx, ∂fx / ∂2fx
    end

    function fdf3(x, p=1)
        fx = f(x, p)
        ∂fx = fp(x)
        ∂2fx = fpp(x)
        ∂3fx = fppp(x)
        return fx, fx / ∂fx, ∂fx / ∂2fx, ∂2fx / ∂3fx
    end

    function fdf4(x, p=1)
        fx = f(x, p)
        ∂fx = fp(x)
        ∂2fx = fpp(x)
        ∂3fx = fppp(x)
        ∂4fx = fpppp(x)
        return fx, fx / ∂fx, ∂fx / ∂2fx, ∂2fx / ∂3fx, ∂3fx / ∂4fx
    end

    @test solve(ZeroProblem((f,), x₀), Roots.LithBoonkkampIJzerman(3, 0)) ≈ α
    @test solve(ZeroProblem(f, x₀), Roots.LithBoonkkampIJzerman(3, 0)) ≈ α
    @test solve(ZeroProblem((f,), x₀), Roots.LithBoonkkampIJzerman(4, 0)) ≈ α
    @test solve(ZeroProblem(f, x₀), Roots.LithBoonkkampIJzerman(4, 0)) ≈ α
    @test solve(ZeroProblem((f,), x₀), Roots.LithBoonkkampIJzerman(3, 0), 1) ≈ α
    @test solve(ZeroProblem(f, x₀), Roots.LithBoonkkampIJzerman(3, 0), 1) ≈ α

    @test solve(ZeroProblem((f, fp), x₀), Roots.LithBoonkkampIJzerman(2, 1)) ≈ α
    @test solve(ZeroProblem((f, fp), x₀), Roots.LithBoonkkampIJzerman(3, 1)) ≈ α
    @test solve(ZeroProblem(fdf1, x₀), Roots.LithBoonkkampIJzerman(2, 1)) ≈ α
    @test solve(ZeroProblem(fdf1, x₀), Roots.LithBoonkkampIJzerman(2, 1), 1) ≈ α

    @test solve(ZeroProblem((f, fp, fpp), x₀), Roots.LithBoonkkampIJzerman(1, 2)) ≈ α
    @test solve(ZeroProblem((f, fp, fpp), x₀), Roots.LithBoonkkampIJzerman(2, 2)) ≈ α
    @test solve(ZeroProblem(fdf2, x₀), Roots.LithBoonkkampIJzerman(1, 2)) ≈ α
    @test solve(ZeroProblem(fdf2, x₀), Roots.LithBoonkkampIJzerman(1, 2), 1) ≈ α

    @test solve(ZeroProblem((f, fp, fpp, fppp), x₀), Roots.LithBoonkkampIJzerman(1, 3)) ≈ α
    @test solve(ZeroProblem((f, fp, fpp, fppp), x₀), Roots.LithBoonkkampIJzerman(2, 3)) ≈ α
    @test solve(ZeroProblem(fdf3, x₀), Roots.LithBoonkkampIJzerman(1, 3)) ≈ α
    @test solve(ZeroProblem(fdf3, x₀), Roots.LithBoonkkampIJzerman(1, 3), 1) ≈ α

    @test solve(ZeroProblem(fdf4, x₀), Roots.LithBoonkkampIJzerman(1, 4)) ≈ α
    @test solve(ZeroProblem(fdf4, x₀), Roots.LithBoonkkampIJzerman(1, 4), 1) ≈ α

    @test solve(
        ZeroProblem((f, fp, fpp, fppp, fpppp), x₀),
        Roots.LithBoonkkampIJzerman(1, 4),
    ) ≈ α
    @test solve(
        ZeroProblem((f, fp, fpp, fppp, fpppp), x̃₀),
        Roots.LithBoonkkampIJzerman(2, 4),
    ) ≈ α # needs closer

    # bracketing
    @test solve(ZeroProblem((f, fp), (1, 2)), Roots.LithBoonkkampIJzermanBracket()) ≈ α
    @test solve(ZeroProblem(fdf1, (1, 2)), Roots.LithBoonkkampIJzermanBracket()) ≈ α
    @test solve(ZeroProblem(fdf1, (1, 2)), Roots.LithBoonkkampIJzermanBracket(), 1) ≈ α
end
