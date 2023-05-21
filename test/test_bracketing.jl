using Roots
using Test
using Printf

## testing bracketing methods
## Originally by John Travers
#
#
## This set of tests is very useful for benchmarking the number of function
## calls, failures, and max errors for the various bracketing methods.
## Table 1 from TOMS748 by Alefeld, Potra, Shi

mutable struct Func
    name::Symbol
    val::Function
    bracket::Function
    params::Vector{Any}
end
function show(io::IO, f::Func)
    @printf io "Func(%s)" f.name
end

## Construct a function object, and check root brackets
macro Func(name)
    @gensym f p b
    esc(quote
        $f = Func($name, val, bracket, params)
        for $p in params
            $b = bracket($p)
            @assert val($p, $b[1]) * $f.val($p, $b[2]) < 0 "Invalid bracket"
        end
        push!(known_functions, $f)
        $f
    end)
end

known_functions = Func[]

## This set of tests is very useful for benchmarking the number of function
## calls, failures, and max errors for the various bracketing methods.

## Table 1 from TOMS748 by Alefeld, Potra, Shi

func1 = let
    val = (_, x) -> sin(x) - x / 2
    bracket(_) = [0.5pi, pi]
    params = [()]
    @Func :func1
end

func2 = let
    val = (n, x) -> -2 * sum([(2i - 5)^2 / (x - i * i)^3 for i in 1:20])
    bracket(n) = [n^2 + 1e-9, (n + 1)^2 - 1e-9]
    params = 1:10
    @Func :func2
end

func3 = let
    val = (p, x) -> p[1] * x * exp(p[2] * x)
    bracket(p) = [-9.0, 31.0]
    params = [(-40.0, -1.0), (-100.0, -2.0), (-200.0, -3.0)]
    @Func :func3
end

func4 = let
    val = (p, x) -> x^p[2] - p[1]
    bracket(p) = p[3]
    params = Tuple{Float64,Float64,Vector{Float64}}[]
    for a_ in [0.2, 1.0], n in 4:2:12
        push!(params, (a_, n, [0.0, 5.0]))
    end
    for n in 8:2:14
        push!(params, (1.0, n, [-0.95, 4.05]))
    end
    @Func :func4
end

func5 = let
    val = (p, x) -> sin(x) - 0.5
    bracket(p) = [0.0, 1.5]
    params = [()]
    @Func :func5
end

func6 = let
    val = (n, x) -> 2x * exp(-n) - 2exp(-n * x) + 1.0
    bracket(n) = [0.0, 1.0]
    params = vcat(1:5, 20:20:100)
    @Func :func6
end

func7 = let
    val = (n, x) -> (1 + (1 - n)^2) * x - (1 - n * x)^2
    bracket(n) = [0.0, 1.0]
    params = [5.0, 10.0, 20.0]
    @Func :func7
end

func8 = let
    val = (n, x) -> x^2 - (1 - x)^n
    bracket(n) = [0.0, 1.0]
    params = [2.0, 5.0, 10.0, 15.0, 20.0]
    @Func :func8
end

func9 = let
    val = (n, x) -> (1 + (1 - n)^4) * x - (1 - n * x)^4
    bracket(n) = [0.0, 1.0]
    params = [1.0, 2.0, 4.0, 5.0, 8.0, 15.0, 20.0]
    @Func :func9
end

func10 = let
    val = (n, x) -> exp(-n * x) * (x - 1) + x^n
    bracket(n) = [0.0, 1.0]
    params = [1, 5, 10, 15, 20]
    @Func :func10
end

func11 = let
    val = (n, x) -> (n * x - 1) / ((n - 1) * x)
    bracket(n) = [0.01, 1.0]
    params = [2, 5, 15, 20]
    @Func :func11
end

func12 = let
    val = (n, x) -> x^(1 / n) - n^(1 / n)
    bracket(n) = [1.0, 100.0]
    params = vcat(2:6, 7:2:33)
    @Func :func12
end

func13 = let
    val = (n, x) -> x == 0.0 ? 0.0 : x / exp(1 / (x * x))
    bracket(n) = [-1.0, 4.0]
    params = [()]
    @Func :func13
end

func14 = let
    val = (n, x) -> x >= 0 ? n / 20 * (x / 1.5 + sin(x) - 1) : -n / 20
    bracket(n) = [-1e4, 0.5pi]
    params = 1:40
    @Func :func14
end

func15 = let
    val = (n, x) -> begin
        if x > 2e-3 / (1 + n)
            exp(1) - 1.859
        elseif x < 0
            -0.859
        else
            exp(0.5e3(n + 1)x) - 1.859
        end
    end
    bracket(n) = [-1e4, 1e-4]
    params = vcat(20:40, 100:100:1000)
    @Func :func15
end

mutable struct MethodResults
    name
    evalcount::Int
    maxresidual::Float64
    failures::Vector{Tuple{Func,Int}}
end
MethodResults() = MethodResults(nothing, 0, 0.0, Tuple{Func,Int}[])
show(io::IO, results::MethodResults) = print(
    io,
    "MethodResults($(results.name), evalcount=$(results.evalcount), numfailures=$(length(results.failures)), maxresidual=$(results.maxresidual))",
)

## Run a method on all known functions.
function run_tests(method; verbose=false, trace=false, name=nothing, abandon=false)
    results = MethodResults()
    results.name = name
    for f in known_functions
        for i in 1:length(f.params)
            p = f.params[i]
            evalcount = 0
            function feval(x)
                evalcount += 1
                result = f.val(p, x)
                trace && @printf "%s[%d]: %s ⇒ %s\n" f i x result
                result
            end
            result, residual = nothing, nothing
            try
                result = method(feval, f.bracket(p))
                isnan(result) && error("NaN")
                residual = f.val(p, result)
                verbose &&
                    @printf "%s[%d] ⇒ %d / %s, residual %.5e\n" f i evalcount result residual
            catch ex
                verbose && @printf "%s[%d] ⇒ FAILED: %s\n" f i ex
                push!(results.failures, (f, i))
                abandon && rethrow(ex)
            end
            results.evalcount += evalcount
            ## Some functions might return non-real values on failures
            if isa(result, AbstractFloat) &&
               isa(residual, AbstractFloat) &&
               isfinite(residual)
                results.maxresidual = max(results.maxresidual, abs(residual))
            end
        end
    end
    results
end

## Run a method on all known functions.
function run_test(f, M; verbose=false, trace=false, name=nothing, abandon=false, kwargs...)
    d = DataFrame(
        i=Int[],
        cnt=Int[],
        steps=Int[],
        Δ=Float64[],
        residual=Float64[],
        result=Float64[],
    )
    for (i, p) in enumerate(f.params)
        evalcount = 0
        function feval(x)
            evalcount += 1
            result = f.val(p, x)
            result
        end
        result, residual = nothing, nothing
        l = Roots.Tracks()
        try
            result = find_zero(feval, f.bracket(p), M; tracks=l, kwargs...)
            isnan(result) && error("NaN")
            residual = f.val(p, result)
        catch ex
            result = NaN
            residual = NaN
        end
        Δ = isempty(l.abₛ) ? NaN : -(l.abₛ[end]...)
        push!(d, (i=i, cnt=evalcount, steps=l.steps, Δ=Δ, residual=residual, result=result))
    end
    d
end

function get_stats(M; kwargs...)
    d = DataFrame(
        fn=Int[],
        i=Int[],
        cnt=Int[],
        steps=Int[],
        Δ=Float64[],
        residual=Float64[],
        result=Float64[],
    )
    for (j, fn) in enumerate(known_functions)
        dⱼ = run_test(fn, M; kwargs...)
        dⱼ.fn .= j
        append!(d, dⱼ)
    end
    d = transform(d, 5:7 => ByRow((Δ, ϵ, a) -> min(abs(a)^(-1) * abs(Δ), abs(ϵ))) => :min)
end

# used to test BracketedHalley
# function rt(f, M)
#     for i in 1:length(f.params)
#         p = f.params[i]
#         fn = x ->f.val(p,x)
#         try
#             x = find_zero((fn, fn', fn''), f.bracket(p), M)#; atol=0.0, rtol=0.0)
#                    @show x, fn(x)
#         catch err
#             @show :err
#         end
#     end
# end

@testset "bracketing methods" begin

    ## Test for failures, ideally all of these would be 0
    ## test for residual, ideally small
    ## test for evaluation counts, ideally not so low for these problems

    ## exact_bracket
    Ms = [
        Roots.Brent(),
        Roots.A42(),
        Roots.AlefeldPotraShi(),
        Roots.Chandrapatla(),
        Roots.ITP(),
        Roots.Ridders(),
        Roots.Bisection(),
    ]
    results = [run_tests((f, b) -> find_zero(f, b, M), name="$M") for M in Ms]
    maxfailures = maximum([length(result.failures) for result in results])
    maxresidual = maximum([result.maxresidual for result in results])
    cnts = [result.evalcount for result in results]
    @test maxfailures == 0
    @test maxresidual <= 5e-13
    @test avg(cnts) <= 4700

    ## False position has larger residuals (and failures until maxsteps is increased)
    Ms = [Roots.FalsePosition(i) for i in 1:12]
    results = [run_tests((f, b) -> find_zero(f, b, M), name="$M") for M in Ms]
    maxfailures = maximum([length(result.failures) for result in results])
    maxresidual = maximum([result.maxresidual for result in results])
    cnts = [result.evalcount for result in results]
    @test maxfailures <= 0
    @test maxresidual <= 1e-5
    @test avg(cnts) <= 3000
end

## Some tests for FalsePosition methods
@testset "FalsePosition" begin
    galadino_probs = [
        (x -> x^3 - 1, [0.5, 1.5]),
        (x -> x^2 * (x^2 / 3 + sqrt(2) * sin(x)) - sqrt(3) / 18, [0.1, 1]),
        (x -> 11x^11 - 1, [0.1, 1]),
        (x -> x^3 + 1, [-1.8, 0]),
        (x -> x^3 - 2x - 5, [2.0, 3]),
        ((x, n=5) -> 2x * exp(-n) + 1 - 2exp(-n * x), [0, 1]),
        ((x, n=10) -> 2x * exp(-n) + 1 - 2exp(-n * x), [0, 1]),
        ((x, n=20) -> 2x * exp(-n) + 1 - 2exp(-n * x), [0, 1]),
        ((x, n=5) -> (1 + (1 - n)^2) * x^2 - (1 - n * x)^2, [0, 1]),
        ((x, n=10) -> (1 + (1 - n)^2) * x^2 - (1 - n * x)^2, [0, 1]),
        ((x, n=20) -> (1 + (1 - n)^2) * x^2 - (1 - n * x)^2, [0, 1]),
        ((x, n=5) -> x^2 - (1 - x)^n, [0, 1]),
        ((x, n=10) -> x^2 - (1 - x)^n, [0, 1]),
        ((x, n=20) -> x^2 - (1 - x)^n, [0, 1]),
        ((x, n=5) -> (1 + (1 - n)^4) * x - (1 - n * x)^4, [0, 1]),
        ((x, n=10) -> (1 + (1 - n)^4) * x - (1 - n * x)^4, [0, 1]),
        ((x, n=20) -> (1 + (1 - n)^4) * x - (1 - n * x)^4, [0, 1]),
        ((x, n=5) -> exp(-n * x) * (x - 1) + x^n, [0, 1]),
        ((x, n=10) -> exp(-n * x) * (x - 1) + x^n, [0, 1]),
        ((x, n=20) -> exp(-n * x) * (x - 1) + x^n, [0, 1]),
        ((x, n=5) -> x^2 + sin(x / n) - 1 / 4, [0, 1]),
        ((x, n=10) -> x^2 + sin(x / n) - 1 / 4, [0, 1]),
        ((x, n=20) -> x^2 + sin(x / n) - 1 / 4, [0, 1]),
    ]

    for (fn_, ab) in galadino_probs
        for M in (FalsePosition(i) for i in 1:12)
            g = Cnt(fn_)
            x0_ = find_zero(g, ab, M)
            @test abs(fn_(x0_)) <= 1e-7
            @test g.cnt <= 50
        end
    end
end

@testset "Bracketing edge cases" begin
    Ms = (Bisection(), Roots.A42(), Roots.AlefeldPotraShi())

    # Endpoints can be infinite
    for M in Ms
        @test find_zero(sign, (-Inf, Inf), M) ≈ 0 atol = 1e-16
    end

    # Function can be infinite for Bisection and Float64
    @test @inferred(find_zero(x -> Inf * sign(x - pi), (-Inf, Inf), Bisection())) ≈ pi

    # finds discontinuities, not necessarily zeros
    f = (x, p=0.0) -> 1 / (x - p) #avoid issue with `0` being identified by `_middle`
    for M in Ms
        @test find_zero(f, (-1, 1), M, p=eps()) ≈ eps() atol = 2eps()
    end

    @test iszero(@inferred(find_zero(f, (-1, 1), Roots.Bisection())))
    # XXX changes with relaxed tolerance (adding non-zero xatol)
    #@test_throws Roots.ConvergenceFailed find_zero(f, (-1, 1), Roots.A42())
    #@test_throws Roots.ConvergenceFailed find_zero(f, (-1, 1), Roots.AlefeldPotraShi())

    # subnormals should still be okay

    α = nextfloat(nextfloat(0.0))
    f = x -> x - α
    for M in (Bisection(),) #Ms XXX NOT A42, AlefeldPotraShi with xatol !==0
        @test find_zero(f, (-1, 1), M) == α
    end

    # with NaN, not Inf
    f = x -> abs(x) / x
    for M in Ms
        @test find_zero(f, (-1, 1), M) ≈ 0 atol = eps()
    end

    # points are not evaluated outside boundary; issue #233
    a, b = -1, 1
    f = x -> abs(x) > 1 ? error("out of bounds") : 1.0 - x
    for M in Ms
        @test find_zero(f, (a, b), M) ≈ 1
    end

    f = x -> abs(x) > 1 ? error("out of bounds") : prevfloat(1.0) - x
    for M in Ms
        @test find_zero(f, (a, b), M) ≈ 1
    end

    # check if fa*fb ≥ 0
    for M in (
        Roots.Bisection(),
        Roots.A42(),
        Roots.AlefeldPotraShi(),
        Roots.Brent(),
        Roots.Ridders(),
        Roots.Chandrapatla(),
        Roots.ITP(),
        Roots.FalsePosition(),
    )
        x = find_zero(x -> sin(x), (0, 1))
        @test iszero(x)
        @test_throws ArgumentError find_zero(x -> sin(x), (2, 3)) # no bracket
    end

    # last bit of accuracy, when close issue #368
    @test find_zero(x -> sinpi(-1 / 40 + x / 40) + 1 - x, (0, 2), A42()) == 1.0

    # sloppy bug using isfinite (#373)
    f = x -> 1 - x / (x - 1)^2
    xleft = 1 + eps(BigFloat)
    xright = 3 * xleft
    x = find_zero(f, (xleft, xright))
    @test abs(f(x)) <= 2eps(BigFloat)

    # simple a42()
    m = run_tests(Roots.a42)
    VERSION >= v"1.6" && @test isempty(m.failures)
    @test m.evalcount <= 3000 # paper says 2884, this has 2877
end
