using Test
using JET

Ts = [Float16, Float32, Float64, BigFloat]

f(x) = x^2 - 2
df(x) = 2x
ddf(x) = 2
F = (f, df, ddf)
a, b = 1,2

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
    Roots.LithBoonkkampIJzerman(3, 0),
    Roots.LithBoonkkampIJzerman(4, 0),
    Roots.Sidi(2),
    Roots.Sidi(3)
]

for T ∈ Ts
    for M ∈ derivative_free
        JET.@test_opt find_zero(f, T(a), M)
        JET.@test_call find_zero(f, T(a), M)
        JET.@test_opt find_zero(f, T(a), M; atol=.0001)
        JET.@test_opt solve(ZeroProblem(f, T(a)), M)
        JET.@test_opt solve(ZeroProblem(f, T(a)), M; atol=0.001)
    end
end

bracketing =  [
    Roots.Brent(),
    Roots.A42(),
    Roots.AlefeldPotraShi(),
    Roots.Chandrapatla(),
    Roots.ITP(),
    Roots.Ridders(),
    Roots.Bisection(),
    Roots.FalsePosition(),
    Roots.ModAB(),
]

for T ∈ Ts
    for M ∈ bracketing
        ## work on these!
        (T == BigFloat && M == Roots.Brent()) && continue
        (T == BigFloat && M == Roots.Chandrapatla()) && continue
        (T == BigFloat && M == Roots.ITP()) && continue
        (T == BigFloat && M == Roots.Ridders()) && continue
        (T == BigFloat && M == Roots.Bisection()) && continue
        (T == BigFloat && M == Roots.ModAB()) && continue

        JET.@test_opt find_zero(f, (T(a), T(b)), M)
        JET.@test_opt find_zero(f, (T(a), T(b)), M; atol=.0001)
        JET.@test_opt solve(ZeroProblem(f, (T(a), T(b))), M)
        JET.@test_opt solve(ZeroProblem(F, (T(a), T(b))), M; atol=0.001)
    end
end


derivative = [
    Roots.Schroder(),
    Roots.Halley(),
    Roots.Newton(),
    Roots.Order1(),
    Roots.Order0(),
    Roots.QuadraticInverse(),
    Roots.SuperHalley(),
    Roots.ChebyshevLike(),
]


for T ∈ Ts
    for M ∈ derivative
        JET.@test_opt find_zero(F, T(a), M)
        JET.@test_opt find_zero(F, T(a), M; atol=.0001)
        JET.@test_opt solve(ZeroProblem(F, T(a)), M)
        JET.@test_opt solve(ZeroProblem(F, T(a)), M; atol=0.001)
    end
end

#@test_throws Test.FallbackTestSetException JET.@test_opt find_zero(sin, 3.0, Order1(); tracks=Roots.Tracks())
