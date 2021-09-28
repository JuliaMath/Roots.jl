# v1.0 only
# Import dependency.
using BenchmarkTools, Roots, Statistics

#Create benchmark group and benchmarks
benchmarks = BenchmarkGroup()

#Put in specific benchmarks

bracket_methods =
    (bisection=Roots.Bisection(), a42=Roots.A42(), aps=Roots.AlefeldPotraShi())
derivative_free_methods = (
    o0=Roots.Order0(),
    o1=Roots.Order1(),
    o1b=Roots.Order1B(),
    o2=Roots.Order2(),
    o2b=Roots.Order2B(),
    o5=Roots.Order5(),
    o8=Roots.Order8(),
    o16=Roots.Order16(),
)

# collection of doable problems
problems = Dict(
    "f1" => (x -> sin(x), 3.0, (3.0, 4.0)),
    "f2" => (x -> x^5 - x - 1, 1.0, (0.5, 5.0)),
    "f3" => (x -> exp(x) - x^4, 7.0, (5.0, 20.0)),
    "f4" => (x -> cos(x) - x / 2, pi / 4, (0.0, pi / 2)),
    "f5" => (x -> x^2 - exp(x) - 3x + 2, -0.5, (-1.0, 1.0)),
    "f6" => (x -> x^2 - exp(x) - 3x + 2, 2.0, (0.0, 3.0)),
    "f7" => (x -> tanh(x) - tan(x), 7.6, (4.0, 8.0)),
    "f8" => (x -> exp(-x^2 + x + 2) - cos(x) + x^3 + 1, -0.5, (-2.0, 1.0)),
    "f9" => (x -> log(x) + sqrt(x) - 5, 7, (7.0, 10.0)),
    "f10" => (x -> log(x) + sqrt(x) - 5, 20, (7.0, 10.0)),
)

function run_bracket(problems, Ms)
    for (nm, prob) in problems
        fn, x0, ab = prob
        for (mnm, M) in zip(fieldnames(typeof(Ms)), Ms)
            solve(ZeroProblem(fn, ab), M)
            #find_zero(fn, ab, M)
        end
    end
end

function run_bracketing(problems, Ms)
    rts = Float64[]
    for (nm, prob) in problems
        fn, x0, ab = prob
        for M in Ms
            rt = solve(ZeroProblem(fn, ab), M)
            #rt = find_zero(fn, ab, M)
            push!(rts, rt)
        end
    end
    rts
end

function run_derivative_free(problems, Ms)
    rts = Float64[]
    for (nm, prob) in problems
        fn, x0, ab = prob
        for M in Ms
            if M == Order0()
                rt = find_zero(fn, x0, M)
            else
                rt = solve(ZeroProblem(fn, ab), M)
            end

            push!(rts, rt)
        end
    end
    rts
end

function run_simple(problems)
    rts = Float64[]
    for (nm, prob) in problems
        fn, x0, ab = prob
        push!(rts, Roots.bisection(fn, ab[1], ab[2]))
        push!(rts, Roots.bisection(fn, ab[1], ab[2], xatol=1e-6))
        push!(rts, Roots.secant_method(fn, x0))
    end
    rts
end

benchmarks = BenchmarkGroup()

benchmarks["bracketing"] = @benchmarkable run_bracketing($problems, $bracket_methods)
benchmarks["derivative_free"] =
    @benchmarkable run_derivative_free($problems, $derivative_free_methods)
benchmarks["simple"] = @benchmarkable run_simple($problems)

for (nm, prob) in problems
    fn, x0, ab = prob
    @assert fn(ab[1]) * fn(ab[2]) < 0

    Ms = bracket_methods
    for (mnm, M) in zip(fieldnames(typeof(Ms)), Ms)
        benchmarks[nm * "-" * string(mnm)] = @benchmarkable find_zero($fn, $ab, $M)
    end

    Ms = derivative_free_methods
    for (mnm, M) in zip(fieldnames(typeof(Ms)), Ms)
        benchmarks[nm * "-" * string(mnm)] = @benchmarkable find_zero($fn, $x0, $M)
    end

    # simple methods
    u, v = ab
    benchmarks[nm * "-bisection"] = @benchmarkable Roots.bisection($fn, $u, $v)
    benchmarks[nm * "-bisection-atol"] =
        @benchmarkable Roots.bisection($fn, $u, $v, xatol=1e-6)
    benchmarks[nm * "-secant"] = @benchmarkable Roots.secant_method($fn, $x0)
end

results = run(benchmarks) # Get results.
results = median(results) # Condense to median.

nm = "benchmarks.json"
fname = joinpath(@__DIR__, nm)
if isinteractive()
    println("""
To save results, manually call in the REPL: BenchmarkTools.save("benchmarks.json", results)
""")
end

#Compare to old results
try
    oldresults = BenchmarkTools.load(fname)[1]
    judge(oldresults, results)
catch err
    error("Couldn't load file- make sure that you've previously saved results.", err.prefix)
end
