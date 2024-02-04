using Roots
using ChainRulesTestUtils
using Zygote
using Test

# issue #325 add frule, rrule
struct 𝐺
    p
end
(g::𝐺)(x) = cos(x) - g.p * x
G₃(p) = find_zero(𝐺(p), (0, pi/2), Bisection())
F₃(p) = find_zero((x,p) -> cos(x) - p*x, (0, pi/2), Bisection(), p)


@testset "Test frule and rrule" begin
    # Type inference tests of `test_frule` and `test_rrule` with the default
    # rule config `ChainRulesTestUtils.TestConfig()` fail due to an issue
    # with ChainRulesTestUtils: https://github.com/JuliaDiff/ChainRulesTestUtils.jl/issues/246

    # single function
    f(x, p) = log(x) - p
    test_frule(solve, ZeroProblem(f, 1), Order1(), 1.0; check_inferred=false)
    test_rrule(solve, ZeroProblem(f, 1), Order1(), 1.0; check_inferred=false)
    if isdefined(Zygote, :ZygoteRuleConfig)
        test_rrule(
            Zygote.ZygoteRuleConfig(),
            solve,
            ZeroProblem(f, 1),
            Order1(),
            1.0;
            check_inferred=false,
        )
    end
    F(p) = find_zero(f, 1, Order1(), p)
    @test first(Zygote.gradient(F, 1)) ≈ exp(1)

    g(x, p) = x^2 - p[1] * x - p[2]
    test_frule(solve, ZeroProblem(g, 1), Order1(), [0.0, 4.0]; check_inferred=false)
    test_rrule(solve, ZeroProblem(g, 1), Order1(), [0.0, 4.0]; check_inferred=false)
    G(p) = find_zero(g, 1, Order1(), p)
    @test first(Zygote.gradient(G, [0, 4])) ≈ [1 / 2, 1 / 4]

    # a tuple of functions
    fx(x, p) = 1 / x
    test_frule(solve, ZeroProblem((f, fx), 1), Roots.Newton(), 1.0; check_inferred=false)
    test_rrule(solve, ZeroProblem((f, fx), 1), Roots.Newton(), 1.0; check_inferred=false)
    if isdefined(Zygote, :ZygoteRuleConfig)
        test_rrule(
            Zygote.ZygoteRuleConfig(),
            solve,
            ZeroProblem((f, fx), 1),
            Roots.Newton(),
            1.0;
            check_inferred=false,
        )
    end
    F2(p) = find_zero((f, fx), 1, Roots.Newton(), p)
    @test first(Zygote.gradient(F2, 1)) ≈ exp(1)

    gx(x, p) = 2x - p[1]
    test_frule(
        solve,
        ZeroProblem((g, gx), 1),
        Roots.Newton(),
        [0.0, 4.0];
        check_inferred=false,
    )
    test_rrule(
        solve,
        ZeroProblem((g, gx), 1),
        Roots.Newton(),
        [0.0, 4.0];
        check_inferred=false,
    )
    G2(p) = find_zero((g, gx), 1, Roots.Newton(), p)
    @test first(Zygote.gradient(G2, [0, 4])) ≈ [1 / 2, 1 / 4]

    # test Functor; issue #408
    x = rand()
    @test first(Zygote.gradient(F₃, x)) ≈ first(Zygote.gradient(G₃, x))
    @test_broken first(Zygote.hessian(F₃, x)) ≈ first(Zygote.hessian(G₃, x))
    # test_frule, test_rrule aren't successful
    #=
    # DimensionMismatch: arrays could not be broadcast to a common size; got a dimension with lengths 3 and 2
    test_frule(
        solve,
        ZeroProblem(𝐺(2), (0.0, pi/2)),
        Roots.Bisection();
        check_inferred=false,
    )
    # MethodError: no method matching keys(::NoTangent)
    test_rrule(
        solve,
        ZeroProblem(𝐺(2), (0.0, pi/2)),
        Roots.Bisection();
        check_inferred=false,
    )
    =#
end
