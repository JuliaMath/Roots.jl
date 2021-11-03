## tests of find_zero interface
using Roots
using Test
using ForwardDiff;
Base.adjoint(f) = x -> ForwardDiff.derivative(f, float(x));

# for a user-defined method
import Roots.Setfield
import Roots.Setfield: @set!
struct Order3_Test <: Roots.AbstractSecant end

## Test the interface
@testset "find_zero interface tests" begin
    meths = [
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
    ]

    ## different types of initial values
    for m in meths
        @test find_zero(sin, 3, m) ≈ pi
        @test find_zero(sin, 3.0, m) ≈ pi
        @test find_zero(sin, big(3), m) ≈ pi
        @test find_zero(sin, big(3.0), m) ≈ pi
    end

    ## defaults for method argument
    @test find_zero(sin, 3.0) ≈ pi    # order0()
    @test @inferred(find_zero(sin, (3, 4))) ≈ π   # Bisection()
    @test @inferred(find_zero(sin, [3, 4])) ≈ π   # Bisection()

    ## test tolerance arguments
    ## xatol, xrtol, atol, rtol, maxevals, maxfneval, strict
    fn, xstar = x -> sin(x) - x + 1, 1.9345632107520243
    x0, M = 20.0, Order2()
    @test find_zero(fn, x0, M) ≈ xstar   # needs 16 iterations, 33 fn evaluations, difference is exact

    # test of maxevals
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M, maxevals=2)
    # test of maxfneval REMOVED
    # @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M, maxfnevals=2)

    # tolerance on f, atol, rtol: f(x) ~ 0
    M = Order2()
    h = 1e-2
    rt = find_zero(fn, x0, M, atol=h, rtol=0.0)
    @test abs(fn(rt)) > h^2 / 100
    rt = find_zero(fn, x0, M, atol=0.0, rtol=h)
    @test abs(fn(rt)) > h^2 / 100

    ## test of tolerances xatol, xrtol with bisection
    a, b = 1.5, 2.0
    h = 1e-6
    M = Roots.Bisection()
    tracks = Roots.Tracks(Float64, Float64)
    @inferred(find_zero(fn, (a, b), M, tracks=tracks, xatol=h, xrtol=0.0))
    u, v = tracks.xs[(end - 1):end]
    @test h >= abs(u - v) >= h / 2

    ## test of strict
    fn, x0 = x -> cos(x) - 1, pi / 4
    @test fn(find_zero(fn, x0, Order5())) <= 1e-8
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, Order5(), strict=true)

    # xn increment needs atol setting for zeros near 0.0 if strict=true
    M = Order1()
    fn = x -> x * exp(x) + nextfloat(0.0)
    @test_throws Roots.ConvergenceFailed find_zero(
        fn,
        1.0,
        M,
        atol=0.0,
        rtol=0.0,
        strict=true,
        xatol=0.0,
    )
    @test abs(find_zero(fn, 1.0, M, atol=0.0, rtol=0.0, strict=true)) <= eps()

    ## test of extreme values for fn, bisection
    c = pi
    fn = x -> Inf * sign(x - c)
    @inferred(find_zero(fn, (-Inf, Inf))) ≈ c

    fn = x -> Inf * x / abs(x) # stop at NaN values
    @inferred(find_zero(fn, (-Inf, Inf))) ≈ 0

    bracketing_meths = (
        Roots.Bisection(),
        Roots.A42(),
        Roots.AlefeldPotraShi(),
        Roots.Brent(),
        Roots.FalsePosition(),
        Roots.FalsePosition(2),
    )

    # test flexbility in interval specification
    for M in bracketing_meths
        @test @inferred(find_zero(sin, (3, 4))) ≈ pi
        @test @inferred(find_zero(sin, [3, 4])) ≈ pi
        @test @inferred(find_zero(sin, 3:4)) ≈ pi
        @test @inferred(find_zero(sin, SomeInterval(3, 4))) ≈ pi
        @test @inferred(find_zero(sin, range(3, stop=4, length=20))) ≈ pi
    end
end

@testset "non simple zeros" begin
    Ms = (
        Roots.Order1B(),
        Roots.Order2B(),
        Roots.Schroder(),
        Roots.Thukral2B(),
        Roots.Thukral3B(),
        Roots.Thukral4B(),
        Roots.Thukral5B(),
    )

    g(x) = exp(x) + x - 2
    f(x) = g(x)^2
    x₀ = 1 / 4

    α = find_zero(g, x₀)
    fs = (f, f', f'', f''', f'''', f''''', f'''''')
    for M in Ms
        @test find_zero(fs, x₀, M) ≈ α atol = 1e-6
    end
end

@testset "find_zero internals" begin

    ##  init_state method
    g1 = x -> x^5 - x - 1
    x0_, xstar_ = (1.0, 2.0), 1.1673039782614187
    M = Roots.A42()
    G1 = Roots.Callable_Function(M, g1)
    state = @inferred(Roots.init_state(M, G1, x0_))
    options = @inferred(Roots.init_options(M, state))
    for M in (Roots.A42(), Roots.Bisection(), Roots.FalsePosition())
        Gₘ = Roots.Callable_Function(M, G1)
        stateₘ = @inferred(Roots.init_state(M, state, Gₘ))
        @test @inferred(solve(M, Gₘ, stateₘ)) ≈ xstar_
    end

    # iterator interface (ZeroProblem, solve; init, solve!)
    meths = [
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
    ]
    g1(x) = x^5 - x - 1
    x0_, xstar_ = 1.16, 1.1673039782614187
    fx = ZeroProblem(g1, x0_)
    for M in meths
        @test solve(fx, M) ≈ xstar_
        P = init(fx, M)
        @test solve!(P) ≈ xstar_
    end

    # solve and parameters
    # should be positional, but named supported for now
    g2 = (x, p) -> cos(x) - x / p
    fx = ZeroProblem(g2, (0, pi / 2))
    @test solve(fx, 2) ≈ @inferred(find_zero(x -> cos(x) - x / 2, (0, pi / 2)))
    @test solve(fx, p=2) ≈ @inferred(find_zero(x -> cos(x) - x / 2, (0, pi / 2)))
    @test @inferred(solve(fx, p=3)) ≈ @inferred(find_zero(x -> cos(x) - x / 3, (0, pi / 2)))
    g3 = (x, p) -> cos(x) + p[1] * x - p[2]
    fx = ZeroProblem(g3, (0, pi / 2))
    @test @inferred(solve(fx, p=[-1 / 10, 1 / 10])) ≈
          @inferred(find_zero(x -> cos(x) - x / 10 - 1 / 10, (0, pi / 2)))

    ## test with early evaluation of bracket
    f = x -> sin(x)
    xs = (3.0, 4.0)
    fxs = f.(xs)
    M = Bisection()
    state = @inferred(Roots.init_state(M, f, xs..., fxs..., m=3.5, fm=f(3.5)))
    @test @inferred(find_zero(M, f, state)) ≈ π

    #     ## hybrid
    g1 = x -> exp(x) - x^4
    x0_, xstar_ = (5.0, 20.0), 8.613169456441398
    M = Roots.Bisection()
    G1 = Roots.Callable_Function(M, g1)
    state = @inferred(Roots.init_state(M, G1, x0_))
    options = @inferred(Roots.init_options(M, state, xatol=1 / 2))
    ZPI = @inferred(init(M, G1, state, options))
    ϕ = iterate(ZPI)
    while ϕ !== nothing
        val, st = ϕ
        state, ctr = st
        ϕ = iterate(ZPI, st)
    end

    N = Roots.Order1() # switch to N
    G2 = Roots.Callable_Function(N, G1)
    stateₙ = Roots.init_state(N, state, G2)
    options = Roots.init_options(N, stateₙ)
    x = solve(N, G2, stateₙ, options)
    @test x ≈ xstar_

    ## test creation of new methods
    ## xn - f/f' - f'' * f^2 / 2(f')^3 = xn - r1 - r1^2/r2 is third order,
    # had to previousely define:
    function Roots.update_state(
        M::Order3_Test,
        f,
        o::Roots.AbstractUnivariateZeroState{T,S},
        options,
        l=Roots.NullTracks(),
    ) where {T,S}
        # xn - f/f' - f'' * f^2 / 2(f')^3 = xn - r1 - r1^2/r2 is third order
        xn_1, xn = o.xn0, o.xn1
        fxn_1, fxn = o.fxn0, o.fxn1

        f_10 = (fxn - fxn_1) / (xn - xn_1)
        xn1::T = xn - fxn / f_10
        fxn1::S = f(xn1)

        f01 = (fxn1 - fxn) / (xn1 - xn)

        if isnan(f_10) || iszero(f_10) || isnan(f01) || iszero(f01)
            return (o, true)
        end

        r1 = fxn1 / f01
        r2 = f01 / ((f01 - f_10) / (xn1 - xn_1))

        wn = xn1 - r1 - r1^2 / r2
        fwn::S = f(wn)

        @set! o.xn0 = xn
        @set! o.xn1 = wn
        @set! o.fxn0 = fxn
        @set! o.fxn1 = fwn

        return (o, false)
    end

    g1 = x -> exp(x) - x^4
    @test find_zero(g1, 8.3, Order3_Test()) ≈ find_zero(g1, 8.3, Order1())

    # test many different calling styles
    f(x) = (sin(x), sin(x) / cos(x)) # x -> (f(x), f(x)/f′(x))
    fs(x) = (sin, cos) # (f, f′)
    x0 = (3, 4)
    g(x, p) = begin
        fx = cos(x) - x / p
        (fx, fx / (-sin(x) - 1 / p))
    end
    x0a = (0.0, pi / 2)
    α₂, α₃ = 1.0298665293222589, 1.1701209500026262
    @test find_zero(f, x0) ≈ π
    @test find_zero(f, first(x0)) ≈ π
    @test find_zero(g, x0a, p=2) ≈ α₂
    @test find_zero(g, first(x0a), p=2) ≈ α₂
    Z = ZeroProblem(f, x0)
    Za = ZeroProblem(g, x0a)
    @test solve(Z) ≈ π
    @test solve(Za, 3) ≈ α₃
    @test solve(Za, p=2) ≈ α₂
    @test solve!(init(Z)) ≈ π
    @test solve!(init(Za, 3)) ≈ α₃
    @test solve!(init(Za, p=3)) ≈ α₃
    Ms = (Roots.Secant(), Roots.Bisection(), Roots.Newton())
    for M in Ms
        @test find_zero(f, x0, M) ≈ π
        @test solve(Z, M) ≈ π
        @test solve!(init(Z, M)) ≈ π
        @test find_zero(g, x0a, M, p=2) ≈ α₂
        @test solve(Za, M, 2) ≈ α₂
        @test solve(Za, M, p=2) ≈ α₂
        @test solve!(init(Za, M, 2)) ≈ α₂
    end
end

@testset "find_zero issue tests" begin

    ## Misc tests
    Ms = [
        Order0(),
        Order1(),
        Roots.Order1B(),
        Order2(),
        Roots.Order2B(),
        Order5(),
        Order8(),
        Order16(),
    ]

    ## issues with starting near a maxima. Some bounce out of it, but
    ## one would expect all to have issues
    fn, xstar = x -> x^3 + 4x^2 - 10, 1.365230013414097
    for M in [Order1(), Roots.Order1B(), Order2(), Roots.Order2B(), Order5()]
        @test_throws Roots.ConvergenceFailed find_zero(fn, -1.0, M)
    end
    for M in [Order0(), Roots.Thukral8(), Roots.Thukral16()]
        @test find_zero(fn, -1.0, M) ≈ xstar
    end

    ## non-simple root
    ## convergence can depend on relaxed convergence checked after an issue
    fn, xstar, x0 = x -> cos(x) - 1, 0.0, 0.1
    for M in Ms
        xn = find_zero(fn, x0, M)
        @test abs(fn(xn)) <= 1e-10
    end
    for M in [Order2(), Order5(), Order8(), Order16()]
        @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M, strict=true)
    end

    ## issue with large steps
    fn, x0 = x -> x^20 - 1, 0.5
    for M in Ms[2:end] # not 0, as it uses bracket
        @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M)
    end

    ## issue with large f''
    fn, x0 = cbrt, 1.0
    for M in [Order1(), Order2(), Order5()]
        @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M)
    end
    ### these stop but only because rtol is used for checking f(xn) ~ 0
    for M in [Roots.Thukral8(), Roots.Thukral16()]
        @test abs(find_zero(fn, x0, M) - 0.0) >= 100
    end

    ## similar (http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html)
    function newton_baffler(x)
        a = 1 / 10
        m, b = 1 / 4, 1 / 8

        if x < -a
            m * x - b
        elseif x > a
            m * x + b
        else
            (m * a + b) / a * (x + a) + (-m * a - b)
        end
    end
    for M in
        (Order0(), Order1(), Roots.Order1B(), Order2(), Roots.Order2B(), Order5(), Order8())
        @test abs(find_zero(newton_baffler, 1.0, M)) <= 1e-15
    end
    for M in (Roots.KumarSinghAkanksha(), Roots.Thukral8(), Roots.Thukral16())
        @test_throws Roots.ConvergenceFailed find_zero(newton_baffler, 1.0, M)
    end

    ## Closed issues ###
    ## issue tests: put in tests to ensure closed issues don't reappear.

    ## issue #94; tolerances not matching documentation
    function test_94(; kwargs...)
        g, T = 1.62850, 14.60000
        α, t1, tf = 0.00347, 40.91375, 131.86573
        y, ya, yf = 0.0, 9000.0, 10000.0
        vy = sqrt(2g * (ya - y))
        θ0, θ1 = atan(α * tf), atan(α * (tf - t1))
        I_sintan(x) = tan(x) / 2cos(x) - atanh(tan(x / 2))
        I_sintan(x, y) = I_sintan(y) - I_sintan(x)
        function lhs(θ)
            tRem = (vy - T / α * (sec(θ1) - sec(θ))) / g
            val = -yf + y + vy * tRem - 0.5g * tRem^2 - T / α^2 * I_sintan(θ, θ1)
            val
        end

        M = Roots.FalsePosition()
        x0 = [atan(α * tf), atan(α * (tf - t1))]
        F = Roots.Callable_Function(M, lhs, nothing) #Roots.DerivativeFree(lhs)
        state = Roots.init_state(M, F, x0)
        options = Roots.init_options(M, state)
        l = Roots.Tracks(state)
        find_zero(M, F, state, options, l)

        @test l.steps <= 45 # 15
    end
    test_94()

    ## Issue with quad_step after truncated M-step PR #140
    @test find_zero(x -> tanh(x) - tan(x), 7.36842, Order0()) ≈ 7.068582745628732

    ## Use tolerance on f, not x with bisectoin
    atol = 0.01
    u = @inferred(find_zero(sin, (3, 4), atol=atol))
    @test atol >= abs(sin(u)) >= atol^2

    ## issue #159 bracket with zeros should be found
    @test @inferred(find_zero(x -> x + 1, (-1, 1))) == -1

    ## issue #178 passinig through method
    @test fzero(sin, 3, 4, Roots.Brent()) ≈ π

    ## issue #188 with A42
    f = let a = 0.18
        x -> x * (1 - x^2) / ((x^2 + a^2) * (1 + a^2 * x^2))
    end
    r = 0.05
    xs = (r + 1e-12, 1.0)
    @test @inferred(find_zero(x -> f(r) - f(x), xs, Roots.A42())) ≈ 0.4715797678171889
end

struct _SampleCallableObject end
(::_SampleCallableObject)(x) = x^5 - x - 1

mutable struct Cnt
    cnt::Int
    f
    Cnt(f) = new(0, f)
end
(f::Cnt)(x) = (f.cnt += 1; f.f(x))

@testset "find_zero with other callable types" begin
    Ms = [
        Order0(),
        Order1(),
        Roots.Order1B(),
        Order2(),
        Roots.Order2B(),
        Order5(),
        Order8(),
        Order16(),
    ]

    for M in Ms
        @test find_zero(_SampleCallableObject(), 1.1, M) ≈ 1.1673039782614187
    end

    for M in Ms
        g = Cnt(x -> x^5 - x - 1)
        @test find_zero(g, 1.1, M) ≈ 1.1673039782614187
        @test g.cnt <= 30
    end
end

@testset "find_bracket test" begin
    fs_xs = ((x -> x^5 - x - 1, (0, 2)), (sin, (3, 4)), (x -> exp(x) - x^4, (5, 20)))

    for (f, x0) in fs_xs
        out = Roots.find_bracket(f, x0)
        @test prod(f.(out.bracket)) <= 0
    end

    ## test size of  bracket
    for M in (Roots.BisectionExact(), Roots.A42(), Roots.AlefeldPotraShi())
        for (fn, b) in (
            (x -> x == 0.0 ? 0.0 : x / exp(1 / (x * x)), (-1.0, 4.0)),
            (x -> exp(-15 * x) * (x - 1) + x^15, (0.0, 1.0)),
            (x -> (-40.0) * x * exp(-x), (-9, 31)),
        )
            l = Roots.find_bracket(fn, b, M)
            @test l.exact || abs(l.bracket[2] - l.bracket[1]) <= eps(l.xstar)
        end
    end

    ## subnormal
    x0 = nextfloat(nextfloat(0.0))
    fn = x -> (x - x0)
    b = (0.0, 1.0)
    m = Roots.__middle(b...)  # 1e-154
    M = Roots.BisectionExact()
    l = Roots.find_bracket(fn, b, M)
    ## Here l.exact=true, but  this checks the bracket size
    ab = abs(l.bracket[2] - l.bracket[1])
    #@test  !(ab) <=  eps(l.xstar))
    @test abs(-(l.bracket...)) <= m

    M = Roots.A42()
    l = Roots.find_bracket(fn, b, M)
    ab = abs(l.bracket[2] - l.bracket[1])
    #@test  !(ab <=  eps(l.xstar))
    @test ab <= m

    M = Roots.AlefeldPotraShi()
    l = Roots.find_bracket(fn, b, M)
    ab = abs(l.bracket[1] - l.bracket[1])
    #@test  !(ab <=  eps(l.xstar))
    @test ab ≤ m

    x0 = nextfloat(0.0)
    fn = x -> x <= 0.0 ? x - x0 : x + x0
    b = (-1.0, 1.0)
    M = Roots.BisectionExact()
    l = Roots.find_bracket(fn, b, M)
    @test !l.exact && abs(-(l.bracket...)) <= maximum(eps.(l.bracket))
end

@testset "function evalutions" begin
    function wrapper(f)
        cnt = 0
        x -> begin
            cnt += 1
            f(x)
        end
    end

    # as of v"1.3.0", no more maxfnevals for stopping, just maxevals
    # this is an alternative
    function fz(f, x0::Number, M; maxfnevals=10, kwargs...)
        F = wrapper(f)
        ZPI = init(ZeroProblem(F, x0), M; kwargs...)
        x = NaN * float(x0)
        ϕ = iterate(ZPI)
        while ϕ !== nothing
            x, st = ϕ
            F.cnt.contents >= maxfnevals && return NaN * float(x0)
            ϕ = iterate(ZPI, st)
        end
        x
    end
    f(x) = x^20 - 1
    x0 = 0.9
    M = Order1()
    @test isnan(fz(f, x0, M))  # takes 19 fn evals, not 10

    # test that for update state, fnevals are correctly counted for simpler
    # methods
    fn = (x) -> sin(x)
    x0 = (3, 4)
    M = Order1()
    state = Roots.init_state(M, Roots.Callable_Function(M, fn), x0)
    options = Roots.init_options(M, state)

    for M in (
        Order1(),
        Order2(),
        Order5(),
        Order8(),
        Order16(),
        Roots.Order1B(),
        Roots.Order2B(),
        Roots.BisectionExact(),
        Roots.Brent(),
        Roots.A42(),
        Roots.AlefeldPotraShi(),
    )

        # test initial count
        g = wrapper(fn)
        G = Roots.Callable_Function(M, g)
        Roots.init_state(M, G, x0)
        @test g.cnt.contents ≤ Roots.initial_fncalls(M)

        # test update state
        g = wrapper(fn)
        stateₘ = Roots.init_state(M, state, Roots.Callable_Function(M, f))
        G = Roots.Callable_Function(M, g)
        l = Roots.Tracks(Float64, Float64)
        Roots.update_state(M, G, stateₘ, options, l)
        @test g.cnt.contents == l.fncalls
    end
end

@testset "_extrema" begin
    @test @inferred(Roots._extrema((π, 0))) === (0.0, Float64(π))
    @test @inferred(Roots._extrema([π, 0])) === (0.0, Float64(π))
    @test_throws ArgumentError Roots._extrema(π)
    @test_throws ArgumentError Roots._extrema((π, π))
    @test_throws ArgumentError Roots._extrema([π, π])
end
