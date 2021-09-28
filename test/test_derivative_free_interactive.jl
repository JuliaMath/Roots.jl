## A set of functions for exploring convergence
## in an interactive manner:

### Benchmarking tests
## We have
##
## * visualize_diagonostics(which): to see summaries of the methods
## over the different functions with which in (:summary, :counts,
## :residuals)
##
## * identify_regression(): to identify regressions in counts,
## residuals, failures, or exact answers, as compared with a past
## diagnostic run
##
## * write_out(): to save results of a diagnostic run. Must be run in the Roots/test directory
##
## * compare_convergence(Methods): to generate a table of convergence orders

using Roots
using Test
using JSON
using Printf

## This uses functions defined here:
include(joinpath(@__DIR__, "test_derivative_free.jl"))

mutable struct Wrapper
    f
    n
end
Wrapper(f) = Wrapper(f, 0)
(F::Wrapper)(x) = (F.n += 1; F.f(x))

struct TallyCard
    cnts::Vector{Int}     # fn evals
    res::Vector{Float64} # rts
    rs::Vector{Float64}   # residual
    exact::Vector{Int}    # exact 1, 0, -1 (for fail)
    ds::Vector{Float64}   # difference from alpha(s)
    rng::Vector{Float64}  # range sampled
end
TallyCard() = TallyCard(Int[], Float64[], Float64[], Int[], Float64[], Float64[])

# test a single function, create tally card
function test_function(method, F)
    tc = TallyCard()
    m, M = Inf, -Inf
    for p in F.params
        f = Wrapper(x -> F.val(p, x))
        result, residual = nothing, nothing
        x0 = F.x0(p)
        m = x0 < m ? x0 : m
        M = x0 > M ? x0 : M
        try
            result = method(f, x0)
            cnt = f.n
            residual = f(result)
            d = minimum(abs.(result .- F.alpha))
            exact =
                (
                    iszero(residual) ||
                    f(result) * f(nextfloat(result)) < 0 ||
                    f(result) * f(prevfloat(result)) < 0
                ) ? 1 : 0

            push!(tc.cnts, cnt)
            push!(tc.res, result)
            push!(tc.ds, d)
            push!(tc.rs, residual)
            push!(tc.exact, exact)
        catch err
            push!(tc.cnts, -1)
            push!(tc.res, NaN)
            push!(tc.ds, NaN)
            push!(tc.rs, NaN)
            push!(tc.exact, -1)
        end
    end
    append!(tc.rng, (m, M))
    tc
end

# vector of vectors to an array
function vvta(vs, T=eltype(vs[1]))
    n = length(vs)
    m = length(vs[1])
    A = zeros(T, n, m)
    for j in 1:n
        A[j, :] = vs[j]
    end
    A
end

function vvta1(vs, T)
    n = length(vs)
    m = length(vs[1])
    A = zeros(T, m, n)
    for j in 1:n
        A[:, j] = vs[j]
    end
    A
end

## Return Dict of arrays
function create_diagonostics(Ms, Fs, nms)
    @assert length(nms) == length(Ms)

    out = Array{Any}(undef, length(Ms), length(Fs))

    for (i, M) in enumerate(Ms)
        for (j, F) in enumerate(Fs)
            out[i, j] = test_function(M, F)
        end
    end

    m, n = length(Ms), length(Fs)

    D = Dict(
        "counts" => Dict(),
        "rts" => Dict(),
        "residuals" => Dict(),
        "exact" => Dict(),
        "delta" => Dict(),
        "failed" => Dict(),
        "nms" => nms,
        "size" => (m, n),
    )

    # Counts
    for j in eachindex(Fs)
        jj = string(j) # for JSON serialization
        cnts = [out[i, j].cnts for i in eachindex(Ms)]
        rts = [out[i, j].res for i in eachindex(Ms)]
        rs = [out[i, j].rs for i in eachindex(Ms)]
        exact = [out[i, j].exact for i in eachindex(Ms)]
        ds = [out[i, j].ds for i in eachindex(Ms)]

        M = length(cnts)
        N = div(M, m)

        D["counts"][jj] = vvta(cnts)
        D["rts"][jj] = vvta(rts)
        D["residuals"][jj] = vvta(rs)
        D["exact"][jj] = vvta(exact)
        D["delta"][jj] = vvta(ds)

        # failed
        fc = D["counts"][jj] .< 0
        fd = D["delta"][jj] .> 1.0
        fr = D["residuals"][jj] .> 1.0
        D["failed"][jj] =
            [fc[i, j] || fd[i, j] || fr[i, j] for i in 1:size(fc)[1], j in 1:size(fc)[2]]
    end

    D
end

function visualize_diagnostics(::Val{:summary}, D)

    ## Some diagnostics
    n = length(D["counts"]) # number of functions
    m = size(D["counts"]["1"])[1]

    for i in 1:m
        fs, cnt, exact, maxresidual = 0, 0, 0, 0.0
        for j in 1:n
            jj = string(j)

            fail = D["failed"][jj][i, :]
            cnts = D["counts"][jj][i, :]
            resids = D["residuals"][jj][i, :]
            exacts = D["exact"][jj][i, :]

            fs += sum(fail)
            cnt += sum(cnts[.!fail])
            exact += sum(exacts .== 1)
            rs = resids[.!fail]
            if !isempty(rs)
                maxresidual = max(maxresidual, maximum(rs))
            end
        end
        nm = D["nms"][i]
        print(rpad(nm, 15))
        print(rpad(cnt, 6))
        print(rpad(exact, 6))
        print(rpad(fs, 5))
        println(maxresidual)
    end
end

function visualize_diagnostics(::Val{:counts}, D, j)

    ## Some diagnostics
    n = length(D["counts"]) # number of functions
    m = size(D["counts"]["1"])[1]

    #    for j in 1:n
    jj = string(j)
    counts = D["counts"][jj]
    fails = D["failed"][jj]
    println(" --- Function $jj ---")
    for i in 1:m
        nm = D["nms"][i]
        print(rpad(nm, 15))
        tot, ntot = 0, 0
        for (fail, cnt) in zip(fails[i, :], counts[i, :])
            if fail
                print(lpad(".", 5))
            else
                tot += cnt
                ntot += 1
                print(lpad(cnt, 5))
            end
        end
        avg = round(tot / ntot, digits=1)
        println("\t$avg")
    end
    #    end
end

function visualize_diagnostics(::Val{:residuals}, D, j)

    ## Some diagnostics
    n = length(D["counts"]) # number of functions
    m = size(D["counts"]["1"])[1]

    jj = string(j)
    resids = D["residuals"][jj]
    exacts = D["exact"][jj]
    fails = D["failed"][jj]
    println(" --- Function $jj ---")
    for i in 1:m
        nm = D["nms"][i]
        print(rpad(nm, 15))
        tot, ntot = 0, 0
        for (fail, res, exact) in zip(fails[i, :], resids[i, :], exacts)
            if fail
                ch = "X"
            elseif res > 1e-3
                ch = "x"
            elseif res > 1e-12
                ch = "~"
            elseif exact == 1
                ch = "∘"
            else
                ch = "."
            end
            print(lpad(ch, 2))
        end
        println("")
    end
end

# which in (:all, :summary, :counts, :residuals)
function visualize_diagnostics(D, which=:summary) # :cnts, ...

    ## Some diagnostics
    n = length(D["counts"]) # number of functions
    m = size(D["counts"]["1"])[1]

    println("Visualize diagnostics: :$which ∈ (:all, :summary, :counts, :residuals)\n")

    # summary
    if which in (:all, :summary)
        println("Method         fncnt exact fail maxresidual")
        # one row summarizing each method
        visualize_diagnostics(Val(:summary), D)
    end

    # Counts
    if which in (:all, :counts)
        for j in 1:n
            visualize_diagnostics(Val(:counts), D, j)
        end
    end

    # residuals
    if which in (:all, :residuals)
        println("Residual summary")
        println(
            "Key:\t`∘` - exact\n\t`.` - res < 1e-12\n\t`~` - 1e-12 < res < 1e-3`\n\t`x` - 1e-3 < res\n\t`X` - failed",
        )

        for j in 1:n
            visualize_diagnostics(Val(:residuals), D, j)
        end
    end
end

## write out/ read in summaries to file
function write_out(fname, D)
    io = open(fname, "w")
    JSON.Writer.print(io, D)
    close(io)
end

# annoyingly need to convert values to proper type after
# JSON serializations
function read_in(fname)
    D = JSON.parsefile(fname)
    E = Dict()
    E["nms"] = string.(D["nms"])
    E["size"] = Tuple(D["size"])

    for nm in ("counts", "exact", "failed")
        E[nm] = Dict()
        for (k, v) in D[nm]
            E[nm][k] = vvta1(v, Int)
        end
    end

    E["residuals"] = Dict()
    for (k, v) in D["residuals"]
        for vi in v
            vi[vi .=== nothing] .= NaN
        end
        E["residuals"][k] = vvta1(v, Float64)
    end

    E
end

# compare D to E
function identify_regressions(Dnew, Dold)
    out = String[]

    Dnew["nms"] == Dold["nms"] || return "Names are different"
    Dnew["size"] == Dold["size"] || return "sizes are different"

    for (k, v) in Dnew["counts"]
        A, A1 = v, Dold["counts"][k]
        sum(A[A .> 0]) <= sum(A1[A1 .> 0]) || push!(out, "counts increased for function $k")
        sum(A[A .> 0]) < sum(A1[A1 .> 0]) &&
            push!(out, "✓ counts decreased for function $k")
    end

    for (k, v) in Dnew["exact"]
        A, A1 = v, Dold["exact"][k]
        sum(A .== 1) >= sum(A1 .== 1) || push!(out, "exact decreased for function $k")
        sum(A .== 1) > sum(A1 .== 1) && push!(out, "✓ exact increased for function $k")
    end

    for (k, v) in Dnew["failed"]
        A, A1 = v, Dold["failed"][k]
        sum(A) <= sum(A1) || push!(out, "failed increased for function $k")
        sum(A) < sum(A1) && push!(out, "✓ failed decreased for function $k")
    end

    for (k, v) in Dnew["residuals"]
        A, A1 = v, Dold["residuals"][k]
        for i in eachindex(A)
            newi, oldi = A[i], A1[i]
            if abs(newi) > 1.1 * abs(oldi)
                push!(out, "residuals increased for function $k")
                break
            end
            if abs(newi) < 0.9 * abs(oldi)
                push!(out, "✓ residuals decreased for function $k")
                break
            end
        end
    end

    return out
end

## Main interface for interactive use
fname = joinpath(@__DIR__, "derivative_free_diagnostics.json")
elide_ascii(x, n=12) = length(x) > n ? x[1:(n - 3)] * "..." * x[(end - 1):end] : x
function create_diagonostics()
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

    Ms = [(f, b) -> find_zero(f, b, M) for M in meths] # F(f,b), name
    nms = elide_ascii.([replace(string(M), r"^Roots." => "") for M in meths])
    Fs = known_functions
    create_diagonostics(Ms, Fs, nms)
end
## write out current diagnostic test
function write_out()
    println("Creating diagonostics to save")
    write_out(fname, create_diagonostics())
end

## visualize state
"""
   visualize_diagnostics(which=:summary)

Show diagnostics summary

`which` is one of  `(:all, :summary, :counts, :residuals)`
"""
visualize_diagnostics(which=:summary) = visualize_diagnostics(create_diagonostics(), which)

## identify regressions from currently saved state
"""
    identify_regressions()

Compare current state to saved state.

Use `write_out` to save a state.
"""
function identify_regressions()
    if !isfile(fname)
        println("No previous diagnostic to compare with")
        return String[] # empty
    end

    Dnew = create_diagonostics()
    Dold = read_in(fname)
    out = identify_regressions(Dnew, Dold)

    out
end

## tests for newton, halley
import ForwardDiff: derivative
D(f, n=1) = n > 1 ? D(D(f), n - 1) : x -> derivative(f, float(x))
function derivative_based_diagonostics()
    Ms = (
        (f, b) -> Roots.find_zero((f, D(f)), b, Roots.Newton()),
        (f, b) -> Roots.find_zero((f, D(f), D(f, 2)), b, Roots.Halley()),
        (f, b) -> Roots.find_zero((f, D(f), D(f, 2)), b, Roots.Schroder()),
        (f, b) -> Roots.find_zero((f, D(f)), b, Roots.Order5()),
    )
    nms = ("Newton", "Halley", "Schroder", "Order5")
    Fs = known_functions
    create_diagonostics(Ms, Fs, nms)
end

## Order of convergence
## assuming e_n = x_n - alpha and
## e_{n+1} = C e_n^q this tries to find q by:
## alpha known: e_{n+1}/e_n = (e_n/e_{n-1})^q,
## alpha unknown: f(x_n) - f(alpha) ~ f'(alpha)*(en), so en ~ Cf(x_n)
function COC(M, f, x0, alpha=missing)
    op = precision(BigFloat)
    setprecision(BigFloat, 8 * 256)
    tracks = Roots.Tracks(BigFloat[], BigFloat[])
    try
        find_zero(f, big(x0), M, tracks=tracks)
    catch err
        setprecision(BigFloat, op)
        rethrow()
    end
    setprecision(BigFloat, op)

    if ismissing(alpha)
        fs = tracks.fs
        [
            Float64(log(abs(fs[k + 2] / fs[k + 1])) / log(abs(fs[k + 1] / fs[k]))) for
            k in 1:(length(fs) - 4)
        ]
    else
        xs = tracks.xs
        es = xs .- alpha
        [Float64(log(abs(es[k + 1])) / log(abs(es[k]))) for k in 1:(length(xs) - 3)]
    end
end

### Traditional tests: start nearby and compare convergence
"""

Compare convergences. For each method and for several test functions, computes
* abs(x1-alpha), abs(x2-alpha), abs(x3-alpha);
* the computational order of convergence rc = log|f(xk+1)/f(xk)|/log|f(xk)/f(xk-1)|
* the computed zero
Example
```
compare_convergence((Order1(), Order2(), Order8()))
##
## Example (errors from Sidi, Unified treatment...)
## Let error be En = xn - alpha
## Secant: E_{n+1} = f[x_n, x_{n-1}, alpha] / f[x_n, x_{n-1}] E_n E_{n-1}
## Newton: E_{n+1} = f[xn, xn, alpha]/f'(xn) E_n^2
## Steffensen: E_{n+1} = f[xn, xn+fxn, alpha] / f[xn, xn + fxn] ⋅ (1 f[xn,alpha] ⋅ En^2
using ForwardDiff; D(f) = x -> ForwardDiff.derivative(f, float(x))
struct MNewton <: Roots.AbstractSecant end
function Roots.update_state(M::MNewton, f, o, opts)
    o.xn0, o.fxn0 = o.xn1, o.fxn1
    o.xn1 -= f(o.xn1) / D(f)(o.xn1)
    o.fxn1 = f(o.xn1)
    nothing
end
compare_convergence((Roots.Secant(), MNewton(), Roots.Steffensen()))

```
"""
function compare_convergence(Ms; F=identity)
    fns = (
        (
            x -> log(x^2 + x + 2) - x + 1,
            big"4.4",
            4.152590736757158274996989004767139785813809448259893154635015805935085336704627,
        ),
        (
            x -> sin(x)^2 - x^2 + 1,
            big"1.4",
            1.404491648215341226035086817786868077176602575918625035145218238569654850906246,
        ),
        (x -> exp(-x^2 + x + 2) - cos(x + 1) + x^3 + 1, -big"0.5", -big"1.0"),
        (
            x -> x^11 + x + 1,
            -big"1",
            -8.4439752879202298306029802319765064181966469706324549911391920468474414245176e-01,
        ),
        (x -> (x - 2) * (x^10 + x + 1) * exp(-x - 1), big"1.9", big"2"),
        (
            x -> sin(3x) + x * cos(x),
            big"1",
            1.197769535216271165938579472950989827411047786536025790115116815210444571657156,
        ),
        (
            x -> exp(-x) - 1 + x / 5,
            big"4.5",
            4.96511423174427630369875913132289394405558498679725097281444614478046398795746,
        ),
        (
            x -> exp(sin(x)) - x + 1,
            big"2.3",
            2.630664147927903633975327052350598568584731954733163386430717083451519883744738,
        ),
        (x -> (exp(x^2 + x - 6) - 1) * (x^2 + 1), big"2.2", big"2"),
        (
            x -> x * exp(x^2) - sin(x)^2 + 3cos(x) + 5,
            -big"1.2",
            -1.207647827130918927009416758356084097760235818949538815205924601763336168539887,
        ),
        (x -> x^8 - x^7 + x^4 - x^3 + x - 1, big"1.2", big"1"),
    )

    # formatting fns
    Log10 = x -> begin # formatting
        if iszero(x)
            return (0, 0)
        end
        n = trunc(Int, log10(x)) - 1
        x / (10.0^n), -n
    end
    elide_ascii(x, n=12) = length(x) > n ? x[1:(n - 2)] * ".." * x[(end - 1):end] : x

    for (i, u) in enumerate(fns)
        fn_, x0, xstar = u
        fn = F(fn_)
        for M in Ms
            tracks = Roots.Tracks(BigFloat, BigFloat)
            a = try
                find_zero(fn, x0, M, tracks=tracks)
            catch err
                NaN * x0
            end
            x1, x2, x3 = tracks.xs[1:3]

            coc = log(abs(fn_(x3) / fn_(x2))) / log(abs(fn_(x2) / fn_(x1)))
            X1, X2, X3 = abs.((x1, x2, x3) .- xstar)

            @printf "%s: %s\t%.3f(-%02i)\t%.3f(-%02i)\t%.3f(-%02i)\t%.2f\t%.8f\n" "f$i" elide_ascii(
                replace(string(M), "Roots." => ""),
                8,
            ) Log10(X1)... Log10(X2)... Log10(X3)... coc a
        end
        println("")
    end
end
