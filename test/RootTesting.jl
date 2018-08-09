using Roots

import Base: show
if VERSION >= v"0.7.0-"
    using Printf
end

mutable struct Func
    name :: Symbol
    val :: Function
    bracket :: Function
    params :: Vector{Any}
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
    val = (_, x) -> sin(x) - x/2
    bracket(_) = [0.5pi, pi]
    params = [()]
    @Func :func1
end

func2 = let
    val = (n, x) -> -2*sum([(2i-5)^2/(x-i*i)^3 for i=1:20])
    bracket(n) = [n^2+1e-9, (n+1)^2-1e-9]
    params = 1:10
    @Func :func2
end

func3 = let
    val = (p, x) -> p[1]*x*exp(p[2]*x)
    bracket(p) = [-9., 31.]
    params = [(-40.,-1.), (-100., -2.), (-200., -3.)]
    @Func :func3
end

func4 = let
    val = (p, x) -> x^p[2] - p[1]
    bracket(p) = p[3]
    params = Tuple{Float64, Float64, Vector{Float64}}[]
    for a_ in [0.2, 1.], n in 4:2:12
        push!(params, (a_, n, [0., 5.]))
    end
    for n in 8:2:14
        push!(params, (1., n, [-0.95, 4.05]))
    end
    @Func :func4
end

func5 = let
    val = (p, x) -> sin(x) - 0.5
    bracket(p) = [0., 1.5]
    params = [()]
    @Func :func5
end

func6 = let
    val = (n, x) -> 2x*exp(-n)-2exp(-n*x)+1.
    bracket(n) = [0., 1.]
    params = vcat(1:5, 20:20:100)
    @Func :func6
end

func7 = let
    val = (n, x) -> (1+(1-n)^2)*x-(1-n*x)^2
    bracket(n)= [0., 1.]
    params = [5., 10., 20.]
    @Func :func7
end

func8 = let
    val = (n, x) -> x^2-(1-x)^n
    bracket(n) = [0., 1.]
    params = [2., 5., 10., 15., 20.]
    @Func :func8
end

func9 = let
    val = (n, x) -> (1+(1-n)^4)*x-(1-n*x)^4
    bracket(n) = [0., 1.]
    params = [1., 2., 4., 5., 8., 15., 20.]
    @Func :func9
end

func10 = let
    val = (n, x) -> exp(-n*x)*(x-1) + x^n
    bracket(n) = [0., 1.]
    params = [1, 5, 10, 15, 20]
    @Func :func10
end

func11 = let
    val = (n, x) -> (n*x-1)/((n-1)*x)
    bracket(n) = [0.01, 1.]
    params = [2, 5, 15, 20]
    @Func :func11
end

func12 = let
    val = (n, x) -> x^(1/n)-n^(1/n)
    bracket(n) = [1., 100.]
    params = vcat(2:6, 7:2:33)
    @Func :func12
end

func13 = let
    val = (n, x) -> x == 0. ? 0. : x/exp(1/(x*x))
    bracket(n) = [-1., 4.]
    params = [()]
    @Func :func13
end

func14 = let
    val = (n, x) -> x >= 0 ? n/20*(x/1.5+sin(x)-1) : -n/20
    bracket(n) = [-1e4, 0.5pi]
    params = 1:40
    @Func :func14
end

func15 = let
    val = (n, x) -> begin
        if x > 2e-3/(1+n)
            exp(1) - 1.859
        elseif x < 0
            -0.859
        else
            exp(0.5e3(n+1)x)-1.859
        end
    end
    bracket(n) = [-1e4, 1e-4]
    params = vcat(20:40, 100:100:1000)
    @Func :func15
end

mutable struct MethodResults
    name
    evalcount :: Int
    maxresidual :: Float64
    failures :: Vector{Tuple{Func, Int}}
end
MethodResults() = MethodResults(nothing, 0, 0., Tuple{Func, Int}[])
show(io::IO, results::MethodResults) =
    print(io, "MethodResults($(results.name), evalcount=$(results.evalcount), numfailures=$(length(results.failures)), maxresidual=$(results.maxresidual))")

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
                verbose && @printf "%s[%d] ⇒ %d / %s, residual %.5e\n" f i evalcount result residual
            catch ex
                verbose && @printf "%s[%d] ⇒ FAILED: %s\n" f i ex
                push!(results.failures, (f, i))
                abandon && rethrow(ex)
            end
            results.evalcount += evalcount
            ## Some functions might return non-real values on failures
            if isa(result, AbstractFloat) && isa(residual, AbstractFloat) && isfinite(residual)
                results.maxresidual = max(results.maxresidual, abs(residual))
            end
        end
    end
    results
end

Ms = vcat([Roots.FalsePosition(i) for i in 1:12], Roots.A42(), Roots.AlefeldPotraShi(), Roots.Brent(), Roots.Bisection())
results = [run_tests((f,b) -> find_zero(f, b, M), name="$M") for M in Ms]
@test all([length(result.failures) <= 10 for result in results])

results = [run_tests((f,b) -> find_zero(f, big.(b), M), name="$M") for M in Ms]
@test all([length(result.failures) == 0 for result in results[13:end]])

          


## run but not in general testing suite
function run_benchmark_tests()

    @printf "%s\n" run_tests((f,b) -> find_zero(f, b, Roots.A42()), name="a42")
    @printf "%s\n" run_tests((f,b) -> find_zero(f, b, Bisection()), name="Bisection")
    @printf "%s\n" run_tests((f,b) -> find_zero(f, b, FalsePosition()), name="FalsePosition")


    for m in [Order0(), Order1(), Order2(), Order5(), Order8(), Order16()]
        @printf "%s\n" run_tests((f, b) -> find_zero(f, mean(b), m), name="$m")
    end

    println("---- using BigFloat ----")

    @printf "%s\n" run_tests((f,b) -> find_zero(f, big.(b), Bisection()), name="a42 (no bisection with Big values)")
    @printf "%s\n" run_tests((f,b) -> find_zero(f, big.(b), FalsePosition()), name="FalsePosition")

    for m in [Order0(), Order1(), Order2(), Order5(), Order8(), Order16()]
        @printf "%s\n" run_tests((f, b) -> find_zero(f, mean(big.(b)), m), name="$m/BigFloat")
    end

    @printf "%s\n" run_tests((f, b) -> Roots.secant_method(f, big.(b)), name="secant_method")
    @printf "%s\n" run_tests((f, b) -> Roots.dfree(f, mean(big.(b))), name="dfree")
#    using ForwardDiff
#    D(f,n=1) = n > 1 ? D(D(f), n-1) : x -> ForwardDiff.derivative(f, float(x))
#    @printf "%s\n" run_tests((f, b) -> find_zero(f, D(f), mean(big.(b)), Order5()), name="Order5/D;BigFloat")
#    @printf "%s\n" run_tests((f, b) -> newton(f, D(f), mean(big.(b))), name="newton/BigFloat")
#    @printf "%s\n" run_tests((f, b) -> halley(f, D(f), D(f,2), mean(big.(b))), name="halley/BigFloat")
end


nothing

