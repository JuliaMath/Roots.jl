module BenchmarkBracketingMethods
# from https://github.com/Proektsoftbg/Numerical/blob/main/Numerical-Julia/root_modab_benchmark.jl

using Test
using Roots

# ---------------------------------------------------------------------------
# Function-call counting wrapper
# ---------------------------------------------------------------------------

mutable struct CountedFunc{F} <: Function
    f::F
    count::Int
end
CountedFunc(f) = CountedFunc(f, 0)
(cf::CountedFunc)(x) = (cf.count += 1; cf.f(x))
reset!(cf::CountedFunc) = (cf.count=0; cf)

# ---------------------------------------------------------------------------
# ModAB solver (translated from ModAB.cs)
# ---------------------------------------------------------------------------

function mod_ab(f, left::Real, right::Real, target::Real=0.0; precision::Float64=1e-14)
    x1, x2 = min(left, right), max(left, right)
    y1 = f(x1) - target
    abs(y1) <= precision && return x1

    y2 = f(x2) - target
    abs(y2) <= precision && return x2

    n_max = -Int(floor(log2(precision) / 2.0)) + 1
    eps1 = precision / 4
    eps = precision * (x2 - x1) / 2.0
    if abs(target) > 1
        eps1 *= target
    end

    side = 0
    ans = x1
    bisection = true
    k = 0.25

    for i in 1:200
        if bisection
            x3 = (x1 + x2) / 2.0
            y3 = f(x3) - target
            ym = (y1 + y2) / 2.0
            if abs(ym - y3) < k * (abs(y3) + abs(ym))
                bisection = false
            end
        else
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if x3 < x1 - eps || x3 > x2 + eps
                return NaN
            end
            y3 = f(x3) - target
        end

        # XXX main difference between Roots implementation is
        # tolerances
        if abs(y3) < eps1 || abs(x3 - ans) < eps
            if x1 > x2
                return side == 1 ? x2 : x1
            end
            return clamp(x3, x1, x2)
        end

        ans = x3
        if sign(y1) == sign(y3)
            if side == 1
                m = 1 - y3 / y1
                if m <= 0
                    y2 /= 2
                else
                    y2 *= m
                end
            elseif !bisection
                side = 1
            end
            x1 = x3
            y1 = y3
        else
            if side == -1
                m = 1 - y3 / y2
                if m <= 0
                    y1 /= 2
                else
                    y1 *= m
                end
            elseif !bisection
                side = -1
            end
            x2 = x3
            y2 = y3
        end
        if i % n_max == 0
            bisection = true
        end
    end
    return ans
end

# ---------------------------------------------------------------------------
# Roots.jl solver wrapper
# ---------------------------------------------------------------------------

function make_roots_solver(method, name::String)
    function solver(f, left::Real, right::Real, target::Real=0.0; precision::Float64=1e-14)
        g = target != 0 ? x -> f(x) - target : f
        a, b = min(left, right), max(left, right)
        try
            return find_zero(
                g,
                (a, b),
                method;
                xatol=precision,
                xrtol=0.0,
                #atol=precision, rtol=0.0,
                maxiters=200,
            )
        catch
            return NaN
        end
    end
    return solver
end

bisect_solver = make_roots_solver(Bisection(), "bisect")
brent_solver  = make_roots_solver(Roots.Brent(), "brent")
ridder_solver = make_roots_solver(Roots.Ridders(), "ridder")
a42_solver    = make_roots_solver(Roots.A42(), "A42")
itp_solver    = make_roots_solver(Roots.ITP(), "ITP")
modab_solve   = make_roots_solver(Roots.ModAB(), "ModAB")
# ---------------------------------------------------------------------------
# Problem definition
# ---------------------------------------------------------------------------

struct Problem
    name::String
    f::Function
    a::Float64
    b::Float64
    value::Float64
end
Problem(name, f, a, b) = Problem(name, f, Float64(a), Float64(b), 0.0)

P(x) = x + 1.11111

# ---------------------------------------------------------------------------
# Test problems
# ---------------------------------------------------------------------------

const problems1 = [
    Problem("f01", x -> x^3 - 1, 0.5, 1.5),
    Problem("f02", x -> x^2 * (x^2 / 3 + sqrt(2) * sin(x)) - sqrt(3) / 18, 0.1, 1),
    Problem("f03", x -> 11 * x^11 - 1, 0.1, 1),
    Problem("f04", x -> x^3 + 1, -1.8, 0),
    Problem("f05", x -> x^3 - 2x - 5, 2, 3),
    Problem("f06", x -> 2x * exp(-5) + 1 - 2exp(-5x), 0, 1),
    Problem("f07", x -> 2x * exp(-10) + 1 - 2exp(-10x), 0, 1),
    Problem("f08", x -> 2x * exp(-20) + 1 - 2exp(-20x), 0, 1),
    Problem("f09", x -> (1 + (1 - 5)^2) * x^2 - (1 - 5x)^2, 0, 1),
    Problem("f10", x -> (1 + (1 - 10)^2) * x^2 - (1 - 10x)^2, 0, 1),
    Problem("f11", x -> (1 + (1 - 20)^2) * x^2 - (1 - 20x)^2, 0, 1),
    Problem("f12", x -> x^2 - (1 - x)^5, 0, 1),
    Problem("f13", x -> x^2 - (1 - x)^10, 0, 1),
    Problem("f14", x -> x^2 - (1 - x)^20, 0, 1),
    Problem("f15", x -> (1 + (1 - 5)^4) * x - (1 - 5x)^4, 0, 1),
    Problem("f16", x -> (1 + (1 - 10)^4) * x - (1 - 10x)^4, 0, 1),
    Problem("f17", x -> (1 + (1 - 20)^4) * x - (1 - 20x)^4, 0, 1),
    Problem("f18", x -> exp(-5x) * (x - 1) + x^5, 0, 1),
    Problem("f19", x -> exp(-10x) * (x - 1) + x^10, 0, 1),
    Problem("f20", x -> exp(-20x) * (x - 1) + x^20, 0, 1),
    Problem("f21", x -> x^2 + sin(x / 5) - 1 / 4, 0, 1),
    Problem("f22", x -> x^2 + sin(x / 10) - 1 / 4, 0, 1),
    Problem("f23", x -> x^2 + sin(x / 20) - 1 / 4, 0, 1),
    Problem("f24", x -> (x + 2) * (x + 1) * (x - 3)^3, 2.6, 4.6),
    Problem("f25", x -> (x - 4)^5 * log(x), 3.6, 5.6),
    Problem("f26", x -> (sin(x) - x / 4)^3, 2, 4),
    Problem(
        "f27",
        x -> (81 - P(x) * (108 - P(x) * (54 - P(x) * (12 - P(x))))) * sign(P(x) - 3),
        1,
        3,
    ),
    Problem("f28", x -> sin((x - 7.143)^3), 7, 8),
    Problem("f29", x -> exp((x - 3)^5) - 1, 2.6, 4.6),
    Problem("f30", x -> exp((x - 3)^5) - exp(x - 1), 4, 5),
    Problem("f31", x -> π - 1 / x, 0.05, 5),
    Problem("f32", x -> 4 - tan(x), 0, 1.5),
    Problem("f33", x -> cos(x) - x^3, 0, 4),
    Problem("f34", x -> cos(x) - x, -11, 9),
    Problem("f35", x -> sqrt(abs(x - 2 / 3)) * (x <= 2 / 3 ? 1 : -1) - 0.1, -11, 9),
    Problem("f36", x -> abs(x - 2 / 3)^0.2 * (x <= 2 / 3 ? 1 : -1), -11, 9),
    Problem("f37", x -> (x - 7 / 9)^3 + (x - 7 / 9) * 1e-3, -11, 9),
    Problem("f38", x -> x <= 1 / 3 ? -0.5 : 0.5, -11, 9),
    Problem("f39", x -> x <= 1 / 3 ? -1e-3 : 1 - 1e-3, -11, 9),
    Problem("f40", x -> x == 0 ? 0.0 : 1 / (x - 2 / 3), -11, 9),
    Problem("f41", x -> 2x * exp(-5) - 2exp(-5x) + 1, 0, 10),
    Problem("f42", x -> (x^2 - x - 6) * (x^2 - 3x + 2), 0, π),
    Problem("f43", x -> x^3, -1, 1.5),
    Problem("f44", x -> x^5, -1, 1.5),
    Problem("f45", x -> x^7, -1, 1.5),
    Problem("f46", x -> (exp(-5x) - x - 0.5) / x^5, 0.09, 0.7),
    Problem("f47", x -> 1 / sqrt(x) - 2log(5e3 * sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f48", x -> 1 / sqrt(x) - 2log(5e7 * sqrt(x)) + 0.8, 0.0005, 0.5),
    Problem("f49", x -> x <= 0 ? (-x^3 - x - 1) : (x^(1 / 3) - x - 1), -1, 1),
    Problem("f50", x -> x^3 - 2x - x + 3, -3, 2),
    Problem("f51", x -> log(x), 0.5, 5),
    Problem("f52", x -> (10 - x) * exp(-10x) - x^10 + 1, 0.5, 8),
    Problem("f53", x -> exp(sin(x)) - x - 1, 1.0, 4),
    Problem("f54", x -> 2sin(x) - 1, 0.1, π / 3),
    Problem("f55", x -> (x - 1) * exp(-x), 0.0, 1.5),
    Problem("f56", x -> (x - 1)^3 - 1, 1.5, 3),
    Problem("f57", x -> exp(x^2 + 7x - 30) - 1, 2.6, 3.5),
    Problem("f58", x -> atan(x) - 1, 1.0, 8),
    Problem("f59", x -> exp(x) - 2x - 1, 0.2, 3),
    Problem("f60", x -> exp(-x) - x - sin(x), 0.0, 2),
    Problem("f61", x -> x^2 - sin(x)^2 - 1, -1, 2),
    Problem("f62", x -> sin(x) - x / 2, π / 2, π),
]

const problems2 = [
    Problem("f63", x -> x * exp(x) - 1, -1, 1),
    Problem("f64", x -> tan(x - 1 / 10), -1, 1),
    Problem("f65", x -> sin(x) + 0.5, -1, 1),
    Problem("f66", x -> 4x^5 + x * x + 1, -1, 1),
    Problem("f67", x -> x + x^10 - 1, -1, 1),
    Problem("f68", x -> π^x - ℯ, -1, 1),
    Problem("f69", x -> log(abs(x - 10 / 9)), -1, 1),
    Problem("f70", x -> 1 / 3 + sign(x) * abs(x)^(1 / 3) + x^3, -1, 1),
    Problem("f71", x -> (x + 2 / 3) / (x + 101 / 100), -1, 1),
    Problem("f72", x -> (x * 1e6 - 1)^3, -1, 1),
    Problem("f73", x -> exp(x) * (x * 1e6 - 1)^3, -1, 1),
    Problem("f74", x -> (x - 1 / 3)^2 * atan(x - 1 / 3), -1, 1),
    Problem("f75", x -> sign(3x - 1) * (1 - sqrt(1 - (3x - 1)^2 / 81)), -1, 1),
    Problem("f76", x -> x > (1 - 1e6) / 1e6 ? (1 + 1e6) / 1e6 : -1.0, -1, 1),
    Problem("f77", x -> x != 1 / 21 ? 1 / (21x - 1) : 0.0, -1, 1),
    Problem("f78", x -> x * x / 4 + ceil(x / 2) - 0.5, -1, 1),
    Problem("f79", x -> ceil(10x - 1) + 0.5, -1, 1),
    Problem("f80", x -> x + sin(x * 1e6) / 10 + 1e-3, -1, 1),
    Problem("f81", x -> x > -1 ? 1 + sin(1 / (x + 1)) : -1.0, -1, 1),
    Problem("f82", x -> 202x - 2floor((2x + 1e-2) / 2e-2) - 0.1, -1, 1),
    Problem("f83", x -> (202x - 2floor((2x + 1e-2) / 2e-2) - 0.1)^3, -1, 1),
]

const problems3 = [
    Problem("f84", x -> (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5) - 0.05, 0.5, 5.5),
    Problem("f85", x -> sin(x) - 0.5x - 0.3, -10.0, 10.0),
    Problem("f86", x -> exp(x) - 1 - x - x * x / 2 - 0.005, -2.0, 2.0),
    Problem("f87", x -> 1 / (x - 0.5) - 2 - 0.05, 0.6, 2.0),
    Problem("f88", x -> log(x) - x + 2 - 0.05, 0.1, 3.0),
    Problem("f89", x -> sin(20x) + 0.1x - 0.1, -5.0, 10.0), # <--- 10.0 not 5.0
    Problem("f90", x -> x^3 - 2x^2 + x - 0.025, -1.0, 2.0),
    Problem("f91", x -> x * sin(1 / x) - 0.1 - 0.01, 0.01, 1.0),
]

const problems_oo = [
    Problem("f92", x -> sign(x - eps(0.0)), -Inf, Inf),
    Problem("f92", x -> sign(x - pi), -Inf, Inf),
    Problem("f93", x -> atan(x - eps(0.0)), -Inf, Inf),
    Problem("f93", x -> atan(x - pi), -Inf, Inf),
]

const all_problems = vcat(problems1, problems2, problems3, problems_oo)

# ---------------------------------------------------------------------------
# Solver table
# ---------------------------------------------------------------------------

const solvers = [
    ("bisect", bisect_solver),
    (" brent", brent_solver),
    ("ridder", ridder_solver),
    ("   A42", a42_solver),
    ("   ITP", itp_solver),
    (" ModAB", modab_solve),
    (" modAB", mod_ab),
]

# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------

function run_benchmark(; verbose=false)
    eps = 1e-14
    col_w = 22

    io = IOBuffer()
    # --- Results ---
    println(io, "Results")
    header =
        lpad("Func", 4) * "; " * join([lpad(name, col_w) for (name, _) in solvers], "; ")
    println(io, header)
    for p in all_problems
        line = lpad(p.name, 4) * "; "
        for (name, solver) in solvers
            cf = CountedFunc(p.f)
            try
                if solver === mod_ab
                    result = solver(cf, p.a, p.b, p.value; precision=eps)
                else
                    result = solver(cf, p.a, p.b, p.value; precision=eps)
                end
                s = if isnan(result)
                    lpad("NaN", col_w)
                else
                    lpad(string(result), col_w)
                end
                line *= s * "; "
            catch e
                line *= lpad("ERR", col_w) * "; "
            end
        end
        println(io, line)
    end
    println(io)

    # --- Function evaluation counts ---
    println(io, "Function evaluations")
    header = lpad("Func", 4) * "; " * join([lpad(name, 6) for (name, _) in solvers], "; ")
    println(io, header)
    total = zeros(Int, length(solvers))
    for p in all_problems
        line = lpad(p.name, 4) * "; "
        for (j, (name, solver)) in enumerate(solvers)
            cf = CountedFunc(p.f)
            try
                if solver === mod_ab
                    solver(cf, p.a, p.b, p.value; precision=eps)
                else
                    solver(cf, p.a, p.b, p.value; precision=eps)
                end
                total[j] += cf.count
                line *= lpad(string(cf.count), 6) * "; "
            catch e
                line *= lpad("ERR", 6) * "; "
            end
        end
        println(io, line)
    end

    # Print totals
    line = lpad("SUM", 4) * "; "
    for t in total
        line *= lpad(string(t), 6) * "; "
    end
    println(io, line)
    println(io)
    if verbose
        println(String(take!(io)))
    end
    return first.(solvers), total
end

@testset "count bracketing steps" begin
    nms, cnts = run_benchmark()

    # bisect ≥ A42 ≥ ITP ≥ ModAB
    @test cnts[1] ≥ cnts[4] ≥ cnts[5] ≥ cnts[6]
end

end
