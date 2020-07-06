## tests of find_zero interface
using Roots
using Test


# for a user-defined method
struct Order3_Test <: Roots.AbstractSecant end

## Test the interface
@testset "find_zero interface tests" begin

    meths = [Order0(),
             Order1(), Roots.Order1B(), Roots.King(),
             Order2(), Roots.Steffensen(), Roots.Order2B(), Roots.Esser(),
             Order5(), Roots.KumarSinghAkanksha(),
             Order8(), Roots.Thukral8(),
             Order16(), Roots.Thukral16()]


    ## different types of initial values
    for m in meths
        @test find_zero(sin, 3, m) ≈ pi
        @test find_zero(sin, 3.0, m) ≈ pi
        @test find_zero(sin, big(3), m) ≈ pi
        @test find_zero(sin, big(3.0), m) ≈ pi
    end


    ## defaults for method argument
    @test find_zero(sin, 3.0) ≈ pi    # order0()
    @test find_zero(sin, (3,4)) ≈ π   # Bisection()
    @test find_zero(sin, [3,4]) ≈ π   # Bisection()



    ## test tolerance arguments
    ## xatol, xrtol, atol, rtol, maxevals, maxfneval, strict
    fn, xstar = x -> sin(x) - x + 1, 1.9345632107520243
    x0, M = 20.0, Order2()
    @test find_zero(fn, x0, M)  ≈ xstar   # needs 16 iterations, 33 fn evaluations, difference is exact


    # test of maxevals, maxfneval
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M, maxevals=2)
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, M, maxfnevals=2)


    # tolerance on f, atol, rtol: f(x) ~ 0
    M = Order2()
    h = 1e-2
    rt = find_zero(fn, x0, M, atol=h, rtol=0.0)
    @test abs(fn(rt)) > h^2 / 100
    rt = find_zero(fn, x0, M, atol=0.0, rtol=h)
    @test abs(fn(rt)) > h^2 / 100


    ## test of tolerances xatol, xrtol with bisection
    a,b = 1.5, 2.0
    h = 1e-6
    M = Roots.Bisection()
    tracks = Roots.Tracks(Float64[], Float64[])
    find_zero(fn, (a,b), M, tracks=tracks, xatol=h, xrtol=0.0)
    u,v = tracks.xs[end-1:end]
    @test h >= abs(u-v) >= h/2

    ## test of strict
    fn, x0 = x -> cos(x) - 1, pi/4
    @test fn(find_zero(fn, x0, Order5())) <= 1e-8
    @test_throws Roots.ConvergenceFailed find_zero(fn, x0, Order5(), strict=true)

    # xn increment needs atol setting for zeros near 0.0 if strict=true
    M = Order1()
    @test_throws Roots.ConvergenceFailed find_zero(x -> x * exp(x), 1.0, M, atol=0.0, rtol=0.0, strict=true, xatol=0.0)
    @test abs(find_zero(x -> x * exp(x), 1.0, M, atol=0.0, rtol=0.0, strict=true)) <= eps()


    ## test of extreme values for fn, bisection
    c = pi
    fn = x -> Inf * sign(x - c)
    find_zero(fn, (-Inf, Inf)) ≈ c

    fn = x -> Inf * x/abs(x) # stop at NaN values
    find_zero(fn, (-Inf, Inf)) ≈ 0
end


@testset "find_zero internals" begin
    ##  init_state! method
    g1(x) = x^5 - x - 1
    x0_, xstar_ = (1.0, 2.0), 1.1673039782614187
    M = Roots.A42()
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state)
    for M in (Roots.A42(), Roots.Bisection(), Roots.FalsePosition())
        Roots.init_state!(state, M, g1, x0_)
        @test find_zero(M, g1, options, state) ≈ xstar_
    end


    ## hybrid
    g1(x) = exp(x) - x^4
    x0_, xstar_ = (5.0, 20.0), 8.613169456441398
    M = Roots.Bisection()
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state, xatol=1/2)
    find_zero(M, g1, options, state)
    M = Roots.Order1()
    Roots.init_state!(state, M, g1, state.m[1])
    Roots.init_options!(options, M)
    @test find_zero(M, g1, options, state) ≈ xstar_


    # test conversion between states/options
    Ms = vcat([Roots.FalsePosition(i) for i in 1:12], Roots.A42(), Roots.AlefeldPotraShi(), Roots.Brent(), Roots.Bisection())
    Ns = [Order1(),
          Roots.Order1B(),
          Order2(),
          Roots.Order2B(),
          Order5(),
          Order8(),
          Order16()]

    g1(x) = exp(x) - x^4
    x0_= (5.0, 20.0)
    M = Ms[1]
    state = Roots.init_state(M, g1, x0_)
    options = Roots.init_options(M, state)

    for M in Ms
        @test Roots.init_state!(state, M, g1, x0_) == nothing
        @test Roots.init_options!(options, M) == nothing
        @test Roots.update_state(M, g1, state, options) == nothing
        @test Roots.assess_convergence(M, state, options) == false # none in just 1 step
    end

    x0_ = 8.0
    xstars = find_zeros(g1, -20.0, 20.0)
    for M in Ns
        @test Roots.init_state!(state, M, g1, x0_) == nothing
        @test Roots.init_options!(options, M) == nothing
        @test Roots.update_state(M, g1, state, options)  == nothing
        @test Roots.assess_convergence(M, state, options) == false  # none in 1 step
    end

    ## test creation of new methods
    ## xn - f/f' - f'' * f^2 / 2(f')^3 = xn - r1 - r1^2/r2 is third order,
    # had to previousely define:
    # struct Order3_Test <: Roots.AbstractSecant end
    function Roots.update_state(M::Order3_Test, f, o::Roots.UnivariateZeroState{T,S}, options) where {T, S}
        # xn - f/f' - f'' * f^2 / 2(f')^3 = xn - r1 - r1^2/r2 is third order
        xn_1, xn = o.xn0, o.xn1
        fxn_1, fxn = o.fxn0, o.fxn1

        f_10 = (fxn - fxn_1) / (xn - xn_1)
        xn1::T = xn - fxn / f_10
        fxn1::S = f(xn1)

        f01 = (fxn1 - fxn) /  (xn1 - xn)

        if isnan(f_10) || iszero(f_10) || isnan(f01) || iszero(f01)
            o.message = "issue with bracket. "
            o.stopped = true
            return
        end

        r1 = fxn1 / f01
        r2 = f01 / ((f01 - f_10) / (xn1 - xn_1))

        wn = xn1 - r1 - r1^2/r2
        fwn::S = f(wn)

        o.xn0, o.xn1 = xn, wn
        o.fxn0, o.fxn1 = fxn, fwn
    end

    g1(x) = exp(x) - x^4
    @test find_zero(g1, 8.3, Order3_Test()) ≈ find_zero(g1, 8.3, Order1())

end

@testset "find_zero issue tests" begin

    ## Misc tests
    Ms = [Order0(), Order1(), Roots.Order1B(), Order2(), Roots.Order2B(),
          Order5(), Order8(), Order16()]

    ## issues with starting near a maxima. Some bounce out of it, but
    ## one would expect all to have issues
    fn, xstar = x -> x^3 + 4x^2 -10,  1.365230013414097
    for M in [Order1(), Roots.Order1B(), Order2(), Roots.Order2B(), Order5()]
        @test_throws Roots.ConvergenceFailed find_zero(fn, -1.0, M)
    end
    for M in [Order0(),  Order8(), Order16()]
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
    for M in [Order8(), Order16()]
        @test abs(find_zero(fn, x0, M) -  0.0) >= 100
    end

    ## similar (http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html)
    function newton_baffler(x)
        a = 1/10
        m, b = 1/4, 1/8

        if x < -a
            m*x - b
        elseif x > a
            m*x + b
        else
            (m*a+b)/a * (x+a) + (-m*a -b)
        end
    end
    for M in (Order0(), Order1(), Roots.Order1B(), Order2(), Roots.Order2B())
        @test abs(find_zero(newton_baffler, 1.0, M)) <= 1e-15
    end
    for M in (Order5(), Order8())
        @test_throws Roots.ConvergenceFailed find_zero(newton_baffler, 1.0, M)
    end


    ## Closed issues ###
    ## issue tests: put in tests to ensure closed issues don't reappear.

    ## issue #94; tolerances not matching documentation
    function test_94(;kwargs...)
        g, T = 1.62850, 14.60000
        α, t1, tf = 0.00347, 40.91375, 131.86573
        y, ya, yf = 0.0, 9000.0, 10000.0
        vy = sqrt(2g*(ya-y))
        θ0, θ1 = atan(α*tf), atan(α*(tf-t1))
        I_sintan(x) = tan(x)/2cos(x) - atanh(tan(x/2))
        I_sintan(x, y) = I_sintan(y) - I_sintan(x)
        function lhs(θ)
        tRem = (vy - T/α*(sec(θ1) - sec(θ))) / g
            val = -yf + y + vy*tRem - 0.5g*tRem^2 - T/α^2*I_sintan(θ, θ1)
            val
        end

        M = Roots.FalsePosition()
        f, x0 = lhs, [atan(α*tf), atan(α*(tf-t1))]
        F = Roots.DerivativeFree(lhs)
        state = Roots.init_state(M, F, x0)
        options = Roots.init_options(M, state)
        find_zero(M, F, options, state)

        @test state.steps <= 45 # 15
    end
    test_94()



    ## Issue with quad_step after truncated M-step PR #140
    @test find_zero(x -> tanh(x) - tan(x), 7.36842 , Order0()) ≈ 7.068582745628732

    ## Use tolerance on f, not x with bisectoin
    atol = 0.01
    u = find_zero(sin, (3, 4), atol=atol)
    @test atol >= abs(sin(u)) >= atol^2

    ## issue #159 bracket with zeros should be found
    @test find_zero(x->x+1,(-1,1)) == -1

    ## issue #178 passinig through method
    @test fzero(sin, 3, 4, Roots.Brent()) ≈ π

    ## issue #188 with A42
    f(x) = x*(1-x^2)/((x^2+a^2)*(1+a^2*x^2))
    a,r = .18, 0.05
    xs = (r + 1e-12, 1.0)
    @test find_zero(x -> f(r)-f(x), xs, Roots.A42()) ≈ 0.4715797678171889
end


struct _SampleCallableObject end
_SampleCallableObject(x) = x^5 - x - 1

mutable struct Cnt
    cnt::Int
    f
    Cnt(f) = new(0, f)
end
(f::Cnt)(x) = (f.cnt += 1; f.f(x))

@testset "find_zero with other callable types" begin

    Ms = [Order0(),
          Order1(), Roots.Order1B(),
          Order2(), Roots.Order2B(),
          Order5(),
          Order8(),
          Order16()]

    for M in Ms
        @test find_zero(_SampleCallableObject, 1.0, M) ≈ 1.1673039782614187
    end

    for M in Ms
        g = Cnt(x -> x^5 - x - 1)
        @test find_zero(g, 1.0, M) ≈ 1.1673039782614187
        @test g.cnt <= 30
    end

end

@testset "find_bracket test" begin
    fs_xs = ((x -> x^5 - x - 1, (0,2)),
             (sin, (3,4)),
             (x -> exp(x) - x^4, (5, 20))
             )

    for (f, x0) in fs_xs
        out = Roots.find_bracket(f, x0)
        @test prod(f.(out.bracket)) <= 0
    end

    ## test size of  bracket
    for M in  (Roots.BisectionExact(),  Roots.A42(), Roots.AlefeldPotraShi())
        for (fn, b) in ((x  -> x == 0. ? 0. : x/exp(1/(x*x)), (-1.0, 4.0)),
                        (x -> exp(-15*x)*(x-1) + x^15, (0., 1.)),
                        (x -> (-40.0)*x*exp(-x),  (-9, 31))
                        )
            l =  Roots.find_bracket(fn, b, M)
            @test  l.exact ||  abs(l.bracket[2] -  l.bracket[1])  <=  eps(l.xstar)
        end
    end

    ## subnormal
    x0 = nextfloat(nextfloat(0.0))
    fn = x  ->  (x-x0)
    b =  (0.0, 1.0)
    m =  Roots.__middle(b...)  # 1e-154
    M = Roots.BisectionExact()
    l = Roots.find_bracket(fn,  b,  M)
    ## Here l.exact=true, but  this checks the bracket size
    @test  !(abs(-(l.bracket...)) <=  eps(l.xstar))
    @test  abs(-(l.bracket...)) <=  m

    M =  Roots.A42()
    l = Roots.find_bracket(fn,  b,  M)
    @test  !(abs(-(l.bracket...)) <=  eps(l.xstar))
    @test  abs(-(l.bracket...)) <= m

    M =  Roots.AlefeldPotraShi()
    l = Roots.find_bracket(fn,  b,  M)
    @test  !(abs(-(l.bracket...)) <=  eps(l.xstar))
    @test  abs(-(l.bracket...)) <= m

    x0 = nextfloat(0.0)
    fn  = x ->  x <=  0.0 ?  x-x0  : x  +  x0
    b = (-1.0, 1.0)
    M = Roots.BisectionExact()
    l = Roots.find_bracket(fn,  b,  M)
    @test !l.exact &&   abs(-(l.bracket...))  <= maximum(eps.(l.bracket))

end
