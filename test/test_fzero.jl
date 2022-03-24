using Test
import Roots.fzero

## Test `fzero` interface to `find_zero`
## test `fzeros` interface for functions

@testset "fzero(s) interface" begin

    ## test interface
    @test fzero(sin, 3) ≈ pi # order0
    @test fzero(sin, (3, 3.1), order=1) ≈ pi # use order if specified
    @test fzero(sin, (3, 4)) ≈ pi  # bracketing
    @test fzero(sin, 3, 4) ≈ pi  # bracketing
    @test fzero(sin, [3, 4]) ≈ pi  # bracketing
    @test_throws ArgumentError fzero(sin, (3, 3.1)) # not a bracket
    @test fzero(sin, cos, 3) ≈ pi # newton

    ## order keyword:
    for o in keys(Roots._method_lookup)
        @test fzero(x -> x^3 - x, 0.9, order=o) ≈ 1.0
    end

    ## bypass order keyword
    for M in
        [Roots.Order0(), Roots.Order1(), Roots.Order1B(), Roots.Order2(), Roots.Order2B()]
        @test fzero(x -> x^3 - x, 0.7, M) ≈ 1.0
    end

    for M in [Roots.Order1(), Roots.Order1B(), Roots.Order2(), Roots.Order2B()]
        N = Roots.Bisection()
        @test fzero(x -> x^3 - x, 0.7, M, N) ≈ 1.0
    end

    ### test tolerances
    fn, xstar, x0, br = x -> x^5 - x - 1, 1.1673039782614187, (1.0, 1.1), [1.0, 2.0]
    @test fzero(fn, x0, order=1) ≈ xstar
    @test !(fzero(fn, x0, order=1, xatol=1e-4, atol=1) ≈ xstar)
    @test_throws Roots.ConvergenceFailed fzero(fn, x0, order=1, maxevals=3)
    @test !(fzero(fn, x0, order=1, maxevals=3, atol=1e-3) ≈ xstar)

    ## Infinities
    ## f(b) = Inf
    f = x -> x + exp(x)
    @test fzero(f, -1e6, 1e6) ≈ -0.5671432904097838
    ## f(a) = Inf
    f = x -> 1 / x - 1
    @test fzero(f, 0, 2) ≈ 1.0

    ## test infinite range
    @test fzero(x -> x, -Inf, Inf) ≈ 0.0

    ##################################################
    ## fzeros function
    ## test interface
    @test fzeros(sin, 3, 4) ≈ [pi]
    @test fzeros(sin, (3, 4)) ≈ [pi]
    @test fzeros(sin, [3, 4]) ≈ [pi]

    @test length(fzeros(x -> exp(x) - x^4, -10, 20)) == 3
    rts = 1:5
    @test all(
        (abs.(fzeros(x -> prod([x - r for r in rts]), 0, 10)) .- collect(1:5)) .<= 1e-15,
    )

    fn = x -> cos(10 * pi * x)
    @test length(fzeros(fn, 0, 1)) == 10

    ### issue with fzeros and roots near 'b'
    @test 0 < maximum(fzeros(x -> sin(x) - 1 / 1000 * x, 0, pi)) < pi
end
