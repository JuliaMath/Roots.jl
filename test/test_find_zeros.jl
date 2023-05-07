# tests for find_zeros

using Roots
using Test

mutable struct CallableFunction
    f
    n
end
(F::CallableFunction)(x) = (F.n += 1; F.f(x))

@testset "find_zeros" begin
    function azero(f, x)
        fx = f(x)
        iszero(fx) && return true
        sign(fx) * sign(f(nextfloat(x))) < 0 && return true
        sign(fx) * sign(f(prevfloat(x))) < 0 && return true
        abs(fx) <= 8eps(x) && return true
        false
    end

    ## easier ones, counting steps
    F = CallableFunction(x -> exp(x) - x^4, 0)
    xrts = find_zeros(F, -5, 20)
    @test length(xrts) == 3
    @test all(azero.((F,), xrts))
    @test F.n <= 1500 #3000

    F = CallableFunction(x -> cos(x) + cos(2x), 0)
    xrts = find_zeros(F, 0, 4pi)
    @test length(xrts) == 6
    @test all(azero.((F,), xrts))
    @test F.n <= 2000 # 5000

    T11(x) = 1024x^11 - 2816x^9 + 2816x^7 - 1232x^5 + 220x^3 - 11x
    U9(x) = 512x^9 - 1024x^7 + 672x^5 - 160x^3 + 10x

    F = CallableFunction(T11, 0)
    xrts = find_zeros(F, -1, 1)
    @test length(xrts) == 11
    @test all(azero.((F,), xrts))
    @test F.n <= 2500 # 10_000

    F = CallableFunction(U9, 0)
    xrts = find_zeros(F, -1, 1)
    @test length(xrts) == 9
    @test all(azero.((F,), xrts))
    @test F.n <= 2500 # 10_000

    W(n) = x -> prod(x - i for i in 1:n)
    Wi(n) = x -> prod((x - i)^i for i in 1:n)

    F = CallableFunction(W(20), 0)
    xrts = find_zeros(F, -1, 21)
    @test length(xrts) == 20
    @test all(azero.((F,), xrts))
    @test F.n <= 4000 #20_000

    F = CallableFunction(Wi(6), 0)
    xrts = find_zeros(F, -1, 7)
    @test length(xrts) == 6
    @test all(azero.((F,), xrts))
    @test F.n <= 10_000

    ## Harder ones
    f1(x) = 2 * exp(0.5 * x) * (sin(5 * x) + sin(101 * x))
    tiger_tail(x) = f1(x) - round(f1(x)) # (-2,1) then filter

    f2(x) = (x - 0.5) * (x - 0.5001) * (x - 4) * (x - 4.05) * (x - 9.3) # (0,10)
    f3(x) = (x - 3)^2 * (x - 4)^2 # (0,5)
    delta = 0.001
    f4(x) = (x - 0.5)^3 * (x - (0.5 + delta)) * (x - 1)
    f5(x) = (x - 0.5)^3 * (x - (0.5 + delta))^3 * (x - 4) * (x - (4 + delta)) * (x - 4.2)^2
    M(n) = x -> prod((x - (0.5 - (1 / 10)^i)) for i in 1:n)
    f6 = M(4)
    f7 = M(5) # too much

    xrts = find_zeros(tiger_tail, -2.0, 1.0)
    xrts = filter(u -> -1 / 4 < tiger_tail(u) < 1 / 4, xrts)
    @test 345 <= length(xrts) <= 345
    @test all(azero.((tiger_tail,), xrts))

    xrts = find_zeros(f2, 0.0, 10.0)
    @test length(xrts) == 5
    @test all(azero.((f2,), xrts))

    xrts = find_zeros(f3, 0.0, 5.0)
    @test length(xrts) == 2
    @test all(azero.((f3,), xrts))

    xrts = find_zeros(f4, 0.0, 1.5)
    @test length(xrts) == 3
    @test all(azero.((f4,), xrts))

    xrts = find_zeros(f5, 0.0, 5.0)
    @test length(xrts) >= 3          # too hard to get 5 w/o luck, as with no_pts=21/k=4
    @test all(azero.((f5,), xrts))

    xrts = find_zeros(f6, 0.0, 10.0)
    @test length(xrts) == 4
    @test all(azero.((f6,), xrts))

    xrts = find_zeros(f7, 0.0, 10.0)  # too sensitive to interval
    @test length(xrts) >= 3           # should be 5
    @test all(azero.((f7,), xrts))

    # Issue #141 solve over [a,b], not (a,b)
    @test length(find_zeros(p -> p * (1 - p), 0, 1)) == 2
    @test length(find_zeros(sin, 0, 5pi)) == 5 + 1

    # test different ways to specify an interval (just need `extrema` defined)
    @test find_zeros(sin, (3, 4)) ≈ [float(pi)]
    @test find_zeros(sin, [3, 4]) ≈ [float(pi)]
    @test find_zeros(sin, 3:4) ≈ [float(pi)]
    @test find_zeros(sin, SomeInterval(3, 4)) ≈ [float(pi)]
    @test find_zeros(sin, range(3, stop=4, length=20)) ≈ [float(pi)]

    # test with constant function
    @test isempty(find_zeros(x -> 4, -10, 10))

    # test with zero function (Issue #339)
    @test_throws DomainError find_zeros(x -> 0, -2, 2)

    # issue #369
    g4(x) = sqrt(abs(x^2 - 1)) / (x * sign(x^2 - 1))
    @test isempty(find_zeros(g4, 1.1, 2))
end

@testset "find_zeros: not Float64 types" begin
    for T in [Float16, Float32, BigFloat]
        rts = find_zeros(x -> cos(x) - x / 10, T(0.0), T(10.0))
        @test eltype(rts) == T
    end
end
