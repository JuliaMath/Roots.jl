using Test
using JET
using Printf


@testset "JET: test package" begin
    JET.test_package(Roots, ignored_modules=(AnyFrameModule(Printf),))
end

Ts = [Float16, Float32, Float64, BigFloat]

f(x) = x^2 - 2
df(x) = 2x
ddf(x) = 2
dddf(x) = 0
F = (f, df, ddf,dddf)
a, b = 1,2

@testset "JET: derivative free" begin
    derivative_free = [
        Order0(),
        Order1(),
        Roots.Order1B(),
        Roots.King(),
        Order2(),
        Roots.Steffensen(),
        Roots.Order2B(),
        Roots.Esser(),
        Order5(),
        Roots.KumarSinghAkanksha(),
        Order8(),
        Roots.Thukral8(),
        Order16(),
        Roots.Thukral16(),
        Roots.LithBoonkkampIJzerman(2, 0),
        Roots.LithBoonkkampIJzerman(3, 0),
        Roots.LithBoonkkampIJzerman(4, 0),
        Roots.LithBoonkkampIJzerman(5, 0),
        Roots.Sidi(2),
        Roots.Sidi(3)
    ]

    for T ∈ Ts
        for M ∈ derivative_free
            JET.@test_opt find_zero(f, T(a), M)
            JET.@test_call find_zero(f, T(a), M)
            JET.@test_opt find_zero(f, T(a), M; atol=.0001)
            JET.@test_call find_zero(f, T(a), M; atol=.0001)
            JET.@test_opt solve(ZeroProblem(f, T(a)), M)
            JET.@test_call solve(ZeroProblem(f, T(a)), M)
            JET.@test_opt solve(ZeroProblem(f, T(a)), M; atol=0.001)
            JET.@test_call solve(ZeroProblem(f, T(a)), M; atol=0.001)
        end
    end
end


@testset "JET: bracketing" begin
    bracketing =  [
        Roots.A42(),
        Roots.AlefeldPotraShi(),
        Roots.Bisection(),
        Roots.Brent(),
        Roots.Chandrapatla(),
        Roots.ITP(),
        Roots.ModAB(),
        Roots.Ridders(),
    ]

    for T ∈ Ts
        for M ∈ bracketing
            if  (T == BigFloat && M == Roots.ModAB()) ||
                (T == BigFloat && M == Roots.ITP())
                # log2 issues with _opt
                JET.@test_call find_zero(f, (T(a), T(b)), M)
                JET.@test_call find_zero(f, (T(a), T(b)), M; atol=.0001)
                JET.@test_call solve(ZeroProblem(f, (T(a), T(b))), M)
                JET.@test_call solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
            else
                JET.@test_opt find_zero(f, (T(a), T(b)), M)
                JET.@test_opt find_zero(f, (T(a), T(b)), M; atol=.0001)
                JET.@test_opt solve(ZeroProblem(f, (T(a), T(b))), M)
                JET.@test_opt solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
                JET.@test_call find_zero(f, (T(a), T(b)), M)
                JET.@test_call find_zero(f, (T(a), T(b)), M; atol=.0001)
                JET.@test_call solve(ZeroProblem(f, (T(a), T(b))), M)
                JET.@test_call solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
            end
        end
    end
end

@testset "JET: FalsePosition" begin
    for T ∈ Ts
        for i ∈ 1:12
            M = FalsePosition(i)
            JET.@test_opt find_zero(f, (T(a), T(b)), M)
            JET.@test_opt find_zero(f, (T(a), T(b)), M; atol=.0001)
            JET.@test_opt solve(ZeroProblem(f, (T(a), T(b))), M)
            JET.@test_opt solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
            JET.@test_call find_zero(f, (T(a), T(b)), M)
            JET.@test_call find_zero(f, (T(a), T(b)), M; atol=.0001)
            JET.@test_call solve(ZeroProblem(f, (T(a), T(b))), M)
            JET.@test_call solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
            end
    end
end

@testset "JET: derivative" begin
    derivative = [
        Roots.Schroder(),
        Roots.Halley(),
        Roots.Newton(),
        Roots.Order1(),
        Roots.Order0(),
        Roots.QuadraticInverse(),
        Roots.SuperHalley(),
        Roots.ChebyshevLike(),
        Roots.Thukral2B(),
        Roots.Thukral3B(),
        Roots.LithBoonkkampIJzerman(1, 1),
        Roots.LithBoonkkampIJzerman(2, 1),
        Roots.LithBoonkkampIJzerman(3, 1),
        Roots.LithBoonkkampIJzerman(1, 2),
        Roots.LithBoonkkampIJzerman(2, 2),
        Roots.LithBoonkkampIJzerman(1, 3),
        Roots.LithBoonkkampIJzerman(2, 3),
    ]

    for T ∈ Ts
        for M ∈ derivative
            JET.@test_opt find_zero(F, T(a), M)
            JET.@test_opt find_zero(F, T(a), M; atol=.0001)
            JET.@test_opt solve(ZeroProblem(F, T(a)), M)
            JET.@test_opt solve(ZeroProblem(F, T(a)), M; atol=0.001)
            JET.@test_call find_zero(F, T(a), M)
            JET.@test_call find_zero(F, T(a), M; atol=.0001)
            JET.@test_call solve(ZeroProblem(F, T(a)), M)
            JET.@test_call solve(ZeroProblem(F, T(a)), M; atol=0.001)

        end
    end
end
