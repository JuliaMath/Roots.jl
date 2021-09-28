using Roots
using Test
using JSON
using Printf

## testing/benchmarking derivative-free methods

mutable struct Func1
    name::Symbol
    val::Function
    x0::Function
    alpha # a value or values
    params::Vector{Any}
end
function Base.show(io::IO, f::Func1)
    @printf io "Func(%s)" f.name
end

#Construct a function object,
macro Func1(name)
    @gensym f p x
    esc(quote
        $f = Func1($name, val, x0, alpha, params)
        push!(known_functions, $f)
        $f
    end)
end

known_functions = Func1[]

## This set of tests is useful for benchmarking the number of function
## calls, failures, and max errors for the various derivative_free methods.

# easy ones
func1 = let
    val = (_, x) -> cos(x) - x / 2
    a, b, N = 0.0, pi / 2, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 1.0298665293222589
    @Func1 :func1
end

func2 = let
    val = (_, x) -> exp(x) - x^4
    a, b, N = 5.0, 20.0, 11
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = find_zeros(x -> val(0, x), -5.0, 20.0)
    @Func1 :func2
end

func3 = let
    val = (_, x) -> x^5 - x - 1
    a, b, N = 0.5, 2.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.1673039782614187
    @Func1 :func3
end

# wider range
func4 = let
    val = (_, x) -> (1 + cos(x)) * (exp(x) - 2)
    a, b, N = 0.0, 1.7, 15
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.6931471805599453
    @Func1 :func4
end

func5 = let
    val = (_, x) -> 2 * x - exp(-x)
    a, b, N = -5.0, 5.0, 15
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.35173371124919584
    @Func1 :func5
end

# once over hump, hard to get back to 0
func6 = let
    val = (_, x) -> x * exp(-x)
    a, b, N = -5.0, 5.0, 15
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.0
    @Func1 :func6
end

# once outside of extrema, points away from zero
func7 = let
    val = (_, x) -> 20.0 * x / (100.0 * x^2 + 1.0)
    a, b, N = -0.2, 0.2, 15
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.0
    @Func1 :func7
end

# should have some failures here

func8 = let
    val = (_, x) -> cos(x) - x / 2
    a, b, N = 0.0, 3pi / 4, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 1.0298665293222589
    @Func1 :func8
end

func9 = let
    val = (_, x) -> tanh(x) - tan(x)
    a, b, N = 5.0, 10.0, 8
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = find_zeros(x -> val(0, x), 5, 12)
    @Func1 :func9
end

func10 = let
    val = (_, x) -> asin(x^2 - 1) - x / 2 + 1
    a, b, N = -0.5, 0.9, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = find_zeros(x -> val(0, x), -0.5, 1.0)
    @Func1 :func10
end

func11 = let
    val = (_, x) -> exp(-x) - cos(x)
    a, b, N = -2.0, 2.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = find_zeros(x -> exp(-x) - cos(x), -2, 2)
    @Func1 :func11
end

func12 = let
    val = (_, x) -> sqrt(x) - 1 / x - 3
    a, b, N = 5, 20, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 9.633595562832696
    @Func1 :func12
end

## multiplicity: http://ir.igsnrr.ac.cn/bitstream/311030/8840/1/%E4%BE%AF%E9%BA%9F%E7%A7%91(SCI)2.pdf
func13 = let
    val = (_, x) -> (x - sqrt(5))^4 / ((x - 1)^2 + 2)
    a, b, N = 3.0, 4.5, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 2.236067977499790
    @Func1 :func13
end

func14 = let
    val = (_, x) -> (sin(x)^2 - 2x + 1)^5
    a, b, N = 3.0, 5.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.71483582544138924
    @Func1 :func14
end

func15 = let
    val = (_, x) -> (8x * exp(-x^2) - 2x - 3)^8
    a, b, N = -2.0, -1.5, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = -1.7903531791589544
    @Func1 :func15
end

# this is 1e6 * (x-1)^7, has many "zeros" near 1
func16 = let
    val = (_, x) -> 1e6 * (x^7 - 7x^6 + 21x^5 - 35x^4 + 35x^3 - 21x^2 + 7x - 1)
    a, b, N = -0.0, 2.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 1.0
    @Func1 :func16
end

func17 = let
    val = (_, x) -> (exp(-x^2 + x + 3) - x + 2)^9
    a, b, N = 2.4, 2.5, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 2.4905398276083051
    @Func1 :func17
end

func18 = let
    val = (_, x) -> (exp(-x) + 2sin(x))^4
    a, b, N = 3.0, 3.4, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 3.1627488709263654
    @Func1 :func18
end

## hard w/o aggressive bracketing
func19 = let
    val = (_, x) -> cbrt(x)
    a, b, N = -3.0, 3.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.0
    @Func1 :func19
end

func20 = let
    val = (_, x) -> x < -1 / 4 ? 3 / 4 * x - 0.3125 : x < 1 / 4 ? 2x : 3 / 4 * x + 0.3125
    a, b, N = -1.0, 1.0, 10
    x0(n) = range(a, stop=b, length=N)[n]
    params = 1:N
    alpha = 0.0
    @Func1 :func20
end

mutable struct MethodResults1
    name
    problems::Int
    evalcount::Int
    maxresidual::Float64
    maxdifference::Float64
    failures::Vector{Tuple{Func1,Int}}
end
MethodResults1() = MethodResults1(nothing, 0, 0, 0.0, Inf, Tuple{Func1,Int}[])
Base.show(io::IO, results::MethodResults1) = print(
    io,
    "MethodResults($(results.name), evalcount=$(results.evalcount), numfailures=$(length(results.failures)), maxresidual=$(results.maxresidual))",
)

## Run a method on all known functions.
mindiff(a, alpha) = minimum([a - i for i in alpha])
function run_tests(method; verbose=false, trace=false, name=nothing, abandon=false)
    results = MethodResults1()
    results.name = name
    results.problems = 0
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
                result = method(feval, f.x0(p))
                md = abs(mindiff(result, f.alpha))
                isnan(result) && error("NaN")
                (md > 1) && error("ran away")

                results.problems += 1
                results.evalcount += evalcount
                residual = f.val(p, result)

                verbose &&
                    @printf "%s[%d] ⇒ %d / %s, residual %.5e\n" f.name i evalcount result residual

                ## Some functions might return non-real values on failures
                if isa(result, AbstractFloat) &&
                   isa(residual, AbstractFloat) &&
                   isfinite(residual)
                    results.maxresidual = max(results.maxresidual, abs(residual))
                    results.maxdifference =
                        max(results.maxdifference, mindiff(result, f.alpha))
                end

                verbose &&
                    abs(residual) > 1 / 10 &&
                    @printf "Large Residual [%s] %d/%d\n" f.name i p

            catch ex
                verbose && @printf "%s[%d] ⇒ FAILED: %s\n" f.name i ex
                push!(results.failures, (f, i))
                abandon && rethrow(ex)
            end
        end
    end
    results
end

avg(x) = sum(x) / length(x)
D(f, h=1e-4) = x -> (f(x + h) - f(x - h)) / (2h)
D2(f, h=1e-4) = x -> (f(x + h) - 2f(x) + f(x - h)) / h^2

if !isinteractive()
    @testset "derivative free methods" begin

        ## Test for failures, ideally all of these would be 0
        ## test for residual, ideally small
        ## test for evaluation counts, ideally not so low for these problems

        ## basic methods
        Ms = [
            Roots.Order1(),
            Roots.Order1B(),
            Roots.Order2(),
            Roots.Order2B(),
            Roots.Order5(),
            Roots.Order8(),
            Roots.Order16(),
        ]
        results = [run_tests((f, b) -> find_zero(f, b, M), name="$M") for M in Ms]

        failures = [length(result.failures) for result in results]
        residuals = [result.maxresidual for result in results]
        cnts = [result.evalcount / result.problems for result in results]

        @test maximum(failures) <= 60
        @test maximum(residuals) <= 5e-5
        @test avg(cnts) <= 40

        ## methods which fall back to bisection when bracket found
        Ms = [Roots.Order0()]
        results = [run_tests((f, b) -> find_zero(f, b, M), name="$M") for M in Ms]

        failures = [length(result.failures) for result in results]
        residuals = [result.maxresidual for result in results]
        cnts = [result.evalcount / result.problems for result in results]

        @test maximum(failures) <= 30
        @test maximum(residuals) <= 1e-5
        @test avg(cnts) <= 40

        ## Newton and Halley
        Fs = [
            (f, b) -> find_zero((f, D(f)), b, Roots.Newton()),
            (f, b) -> find_zero((f, D(f), D2(f)), b, Roots.Halley()),
            (f, b) -> find_zero((f, D(f), D2(f)), b, Roots.Schroder()),
        ]
        results = [run_tests(F) for F in Fs]

        failures = [length(result.failures) for result in results]
        residuals = [result.maxresidual for result in results]
        cnts = [result.evalcount / result.problems for result in results]

        @test maximum(failures) <= 70
        @test maximum(residuals) <= 5e-4
        @test avg(cnts) <= 50
    end

    @testset "derivative free, non Float64" begin
        Ms = [
            Roots.Order0(),
            Roots.Order1(),
            Roots.Order1B(),
            Roots.Order2(),
            Roots.Order2B(),
            Roots.Order5(),
            Roots.Order8(),
            Roots.Order16(),
        ]
        Ts = [Float16, Float32, BigFloat]

        fn = x -> x^3 - 2x - 5
        alpha =
            2.094551481542326591482386540579302963857306105628239180304128529045312189983499
        x0 = 2.0

        for T in Ts
            for M in Ms
                xstar = find_zero(fn, T(x0), M)
                @test xstar ≈ T(alpha) atol = max(sqrt(eps(T)), eps())
                @test isa(xstar, T)
            end
        end

        for T in Ts
            xstar = Roots.find_zero((fn, D(fn, sqrt(eps(T)))), T(x0), Roots.Newton())
            @test xstar ≈ T(alpha) atol = max(sqrt(eps(T)), eps())
            @test isa(xstar, T)

            xstar = Roots.find_zero(
                (fn, D(fn, sqrt(eps(T))), D2(fn, sqrt(eps(T)))),
                T(x0),
                Roots.Halley(),
            )
            @test xstar ≈ T(alpha) atol = max(sqrt(eps(T)), eps())
            @test isa(xstar, T)

            xstar = Roots.find_zero(
                (fn, D(fn, sqrt(eps(T))), D2(fn, sqrt(eps(T)))),
                T(x0),
                Roots.Schroder(),
            )
            @test xstar ≈ T(alpha) atol = max(sqrt(eps(T)), eps())
            @test isa(xstar, T)
        end
    end
end
