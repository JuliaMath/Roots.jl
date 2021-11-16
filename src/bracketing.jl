###

const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""

abstract type AbstractBisection <: AbstractBracketing end
fn_argout(::AbstractBracketing) = 1

"""

    Bisection()
    Roots.BisectionExact()

If possible, will use the bisection method over `Float64` values. The
bisection method starts with a bracketing interval `[a,b]` and splits
it into two intervals `[a,c]` and `[c,b]`, If `c` is not a zero, then
one of these two will be a bracketing interval and the process
continues. The computation of `c` is done by `_middle`, which
reinterprets floating point values as unsigned integers and splits
there. It was contributed  by  [Jason Merrill](https://gist.github.com/jwmerrill/9012954).
This method avoids floating point issues and when the
tolerances are set to zero (the default) guarantees a "best" solution
(one where a zero is found or the bracketing interval is of the type
`[a, nextfloat(a)]`).

When tolerances are given, this algorithm terminates when the interval
length is less than or equal to the tolerance
`max(xtol, max(abs(a), abs(b)) * .xrtol)`.


When a zero tolerance is given and the values are not `Float64`
values, this will call the [`A42`](@ref) method.


"""
struct Bisection <: AbstractBisection end  # either solvable or A42
struct BisectionExact <: AbstractBisection end

## tracks for bisection, different, we show bracketing interval
function log_step(l::Tracks, M::AbstractBracketing, state)
    push!(l.xs, state.xn0)
    push!(l.xs, state.xn1) # we store [ai,bi, ai+1, bi+1, ...]
    log_steps(l)
end
function show_tracks(io::IO, l::Tracks, M::AbstractBracketing)
    xs = l.xs
    n = length(xs)
    for (i, j) in enumerate(1:2:(n - 1))
        println(
            io,
            @sprintf(
                "(%s, %s) = (% 18.16f, % 18.16f)",
                "a_$(i-1)",
                "b_$(i-1)",
                xs[j],
                xs[j + 1]
            )
        )
    end
    println(io, "")
end

## helper function
function adjust_bracket(x0)
    u, v = map(float, _extrema(x0))
    if u > v
        u, v = v, u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end

function init_state(M::AbstractBracketing, F::Callable_Function, x)
    x₀′, x₁′ = map(float, _extrema(x))
    x₀, x₁ = adjust_bracket((x₀′, x₁′))
    fx₀, fx₁ = F(x₀), F(x₁)
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function init_state(::AbstractBracketing, F, x₀, x₁, fx₀, fx₁; m=_middle(x₀, x₁), fm=F(m))
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)

    assert_bracket(fx₀, fx₁)
    if x₀ > x₁
        x₀, x₁, fx₀, fx₁ = x₁, x₀, fx₁, fx₀
    end

    #    xₘ = Roots._middle(x₀, x₁) # for possibly mixed sign x1, x2
    if sign(fm) * fx₀ < 0
        UnivariateZeroState(m, x₀, fm, fx₀)
    else
        UnivariateZeroState(x₁, m, fx₁, fm)
    end
end

initial_fncalls(::Roots.AbstractBracketing) = 3

# for Bisection, the defaults are zero tolerances and strict=true
"""
    default_tolerances(M::AbstractBisection, [T], [S])


For `Bisection` (or `BisectionExact`), when the `x` values are of type `Float64`, `Float32`,
or `Float16`, the default tolerances are zero and there is no limit on
the number of iterations. In this case, the
algorithm is guaranteed to converge to an exact zero, or a point where
the function changes sign at one of the answer's adjacent floating
point values.

For other types,  the [`Roots.A42`](@ref) method (with its tolerances) is used.

"""
function default_tolerances(::AbstractBisection, ::Type{T}, ::Type{S}) where {T,S}
    xatol = zero(T)
    xrtol = zero(one(T))
    atol = zero(float(one(S))) * oneunit(S)
    rtol = zero(float(one(S))) * one(S)
    maxevals = typemax(Int)
    maxfnevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

# find middle of (a,b) with convention that
# * if a, b finite, they are made non-finite
# if a,b of different signs, middle is 0
# middle falls back to a/2 + b/2, but
# for Float64 values, middle is over the
# reinterpreted unsigned integer.
function _middle(x, y)
    a = isinf(x) ? nextfloat(x) : x
    b = isinf(y) ? prevfloat(y) : y
    if sign(a) * sign(b) < 0
        return zero(a)
    else
        __middle(a, b)
    end
end

const FloatNN = Union{Float64,Float32,Float16}

## find middle assuming a,b same sign, finite
## Alternative "mean" definition that operates on the binary representation
## of a float. Using this definition, bisection will never take more than
## 64 steps (over Float64)
__middle(x::Float64, y::Float64) = __middle(Float64, UInt64, x, y)
__middle(x::Float32, y::Float32) = __middle(Float32, UInt32, x, y)
__middle(x::Float16, y::Float16) = __middle(Float16, UInt16, x, y)
## fallback for non FloatNN number types
__middle(x::Number, y::Number) = 0.5 * x + 0.5 * y

function __middle(T, S, x, y)
    # Use the usual float rules for combining non-finite numbers
    # do division over unsigned integers with bit shift
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    sign(x + y) * reinterpret(T, mid)
end

## the method converges,
## as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(::Bisection, state::AbstractUnivariateZeroState, options)
    flag, converged = assess_convergence(BisectionExact(), state, options)
    converged && return (flag, converged)

    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    if isapprox(a, b, atol=options.xabstol, rtol=options.xreltol)
        return (:x_converged, true)
    end

    ftol = max(options.abstol, max(abs(a), abs(b)) * options.reltol)
    if min(abs(fa), abs(fb)) ≤ ftol
        return (:f_converged, true)
    end

    return :not_converged, false
end

# for exact convergence, we can skip some steps
function assess_convergence(::BisectionExact, state::UnivariateZeroState, options)
    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    (iszero(fa) || isnan(fa) || iszero(fb) || isnan(fb)) && return (:f_converged, true)

    nextfloat(a) == b && return (:x_converged, true)

    return (:not_converged, false)
end

# assumes stopped = :x_converged
function decide_convergence(
    ::AbstractBracketing,
    F,
    state::AbstractUnivariateZeroState,
    options,
    val,
)
    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    isnan(fa) && return a
    isnan(fb) && return b

    if abs(fa) < abs(fb)
        return a
    else
        return b
    end
end

##################################################
function update_state(M::AbstractBisection, F, o, options, l=NullTracks())
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1

    c = _middle(a, b)
    fc = F(c)
    incfn(l)

    if sign(fa) * sign(fc) < 0
        b, fb = c, fc
    else
        a, fa = c, fc
    end

    @set! o.xn0 = a
    @set! o.xn1 = b
    @set! o.fxn0 = fa
    @set! o.fxn1 = fb

    return o, false
end

##################################################

## Bisection has special cases
## for zero tolerance, we have either BisectionExact or A42 methods
## for non-zero tolerances, we have use thegeneral Bisection method
function find_zero(
    fs,
    x0,
    method::Bisection;
    p=nothing,
    tracks=NullTracks(),
    verbose=false,
    kwargs...,
)
    _options = init_options(Bisection(), Float64, Float64; kwargs...)
    iszero_tol =
        iszero(_options.xabstol) &&
        iszero(_options.xreltol) &&
        iszero(_options.abstol) &&
        iszero(_options.reltol)

    _method = if iszero_tol
        float(first(_extrema(x0))) isa FloatNN ? BisectionExact() : A42()
    else
        method
    end

    return solve(ZeroProblem(fs, x0), _method, p; verbose=verbose, tracks=tracks, kwargs...)
end

function assess_convergence(
    method::AbstractBracketing,
    state::AbstractUnivariateZeroState,
    options,
)
    a, b, fa, fb = state.xn0, state.xn1, state.fxn0, state.fxn1

    if isnan(a) || isnan(b)
        return (:nan, true)
    end

    if isnan(fa) || isnan(fb)
        return (:nan, true)
    end

    M = maximum(abs, (a, b))
    δₓ = maximum(
        promote(options.xabstol, M * options.xreltol, sign(options.xreltol) * eps(M)),
    )

    if abs(b - a) <= 2δₓ
        return (:x_converged, true)
    end

    # check f
    u, fu = choose_smallest(a, b, fa, fb)
    δ = maximum(promote(options.abstol, M * options.reltol * (oneunit(fu) / oneunit(u))))
    if abs(fu) <= δ
        iszero(fu) && return (:exact_zero, true)
        return (:f_converged, true)
    end

    return (:not_converged, false)
end

## convergence is much different here
function check_zero(::AbstractBracketing, state, c, fc)
    isnan(c) && return true
    isinf(c) && return true
    iszero(fc) && return true
    return false
end

###################################################
#
## Alefeld, Potra, Shi have two algorithms belosw, one is most efficient, but
## slightly slower than other.
abstract type AbstractAlefeldPotraShi <: AbstractBracketing end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided interval `[a, b]`, without requiring derivatives. It is based
on algorithm 4.2 described in: G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995), DOI: [10.1145/210089.210111](https://doi.org/10.1145/210089.210111).
The asymptotic efficiency index, ``q^{1/k}``, is ``(2 + 7^{1/2})^{1/3} = 1.6686...``.
Originally by John Travers.

"""
struct A42 <: AbstractAlefeldPotraShi end

## put in utils?
@inline isbracket(fa, fb) = sign(fa) * sign(fb) < 0
assert_bracket(fx0, fx1) = isbracket(fx0, fx1) || throw(ArgumentError(bracketing_error))

# f[a, b]
@inline f_ab(a, b, fa, fb) = (fb - fa) / (b - a)

# f[a,b,d]
@inline function f_abd(a, b, d, fa, fb, fd)
    fab, fbd = f_ab(a, b, fa, fb), f_ab(b, d, fb, fd)
    (fbd - fab) / (d - a)
end

# assume fc != 0
## return a1,b1,d with a < a1 <  < b1 < b, d not there
@inline function bracket(a, b, c, fa, fb, fc)
    if isbracket(fa, fc)
        # switch b,c
        return (a, c, b, fa, fc, fb)
    else
        # switch a,c
        return (c, b, a, fc, fb, fa)
    end
end

# Cubic if possible, if not, quadratic(3)
function take_a42_step(a, b, d, ee, fa, fb, fd, fe, k, delta=zero(a))
    # if r is NaN or Inf we move on by condition. Faster than checking ahead of time for
    # distinctness
    r = ipzero(a, b, d, ee, fa, fb, fd, fe, delta) # let error and see difference in allcoation?
    (a + 2delta < r < b - 2delta) && return r
    r = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
end

function ipzero(a, b, c, d, fa, fb, fc, fd, delta=zero(a))
    Q11 = (c - d) * fc / (fd - fc)
    Q21 = (b - c) * fb / (fc - fb)
    Q31 = (a - b) * fa / (fb - fa)
    D21 = (b - c) * fc / (fc - fb)
    D31 = (a - b) * fb / (fb - fa)
    Q22 = (D21 - Q11) * fb / (fd - fb)
    Q32 = (D31 - Q21) * fa / (fc - fa)
    D32 = (D31 - Q21) * fc / (fc - fa)
    Q33 = (D32 - Q22) * fa / (fd - fa)
    c = a + (Q31 + Q32 + Q33)

    (a + 2delta < c < b - 2delta) && return c

    newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
end

# return c in (a+delta, b-delta)
# adds part of `bracket` from paper with `delta`
function newton_quadratic(a, b, d, fa, fb, fd, k::Int, delta=zero(a))
    A = f_abd(a, b, d, fa, fb, fd)
    r = isbracket(A, fa) ? b : a

    # use quadratic step; if that fails, use secant step; if that fails, bisection
    if !(isnan(A) || isinf(A)) || !iszero(A)
        B = f_ab(a, b, fa, fb)

        for i in 1:k
            Pr = fa + B * (r - a) + A * (r - a) * (r - b)
            Prp = (B + A * (2r - a - b))
            r -= Pr / Prp
        end
        if a + 2delta < r < b - 2delta
            return r
        end
    end
    # try secant step
    r = secant_step(a, b, fa, fb)

    if a + 2delta < r < b - 2delta
        return r
    end

    return _middle(a, b)
end

# for A42, the defaults are reltol=eps(), atol=0; 45 evals and strict=true
# this *basically* follows the tol in the paper (2|u|*rtol + atol)
"""
    default_tolerances(::AbstractAlefeldPotraShi, T, S)

The default tolerances for Alefeld, Potra, and Shi methods are
`xatol=zero(T)`, `xrtol=eps(T)/2`, `atol= zero(S), and rtol=zero(S)`, with
appropriate units; `maxevals=45`, `maxfnevals = Inf`; and `strict=true`.

"""
default_tolerances(M::AbstractAlefeldPotraShi) = default_tolerances(M, Float64, Float64)
function default_tolerances(::AbstractAlefeldPotraShi, ::Type{T}, ::Type{S}) where {T,S}
    xatol = zero(T)
    xrtol = eps(one(T)) / 2
    atol = zero(float(one(S))) * oneunit(S)
    rtol = zero(float(one(S))) * one(S)
    maxevals = 45
    maxfnevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

## initial step, needs to log a,b,d
function log_step(l::Tracks, M::AbstractAlefeldPotraShi, state, ::Any)
    a, b, c = state.xn0, state.xn1, state.d
    append!(l.xs, extrema((a, b, c)))
    push!(l.xs, a)
    push!(l.xs, b) # we store [ai,bi, ai+1, bi+1, ...] for brackecting methods
    log_steps(l, 1)
end

struct A42State{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    d::T
    ee::T
    fxn1::S
    fxn0::S
    fd::S
    fee::S
end

function init_state(::A42, F, x₀, x₁, fx₀, fx₁; c=_middle(x₀, x₁), fc=F(c))
    a, b, fa, fb = x₀, x₁, fx₀, fx₁
    if a > b
        a, b, fa, fb = b, a, fb, fa
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
    ee, fe = NaN * d, fd # use NaN for initial

    A42State(b, a, d, ee, fb, fa, fd, fe)
end
initial_fncalls(::A42) = 3

# helper: set state xn1, fxn1
function _set_state(state, x₁, fx₁)
    @set! state.xn1 = x₁
    @set! state.fxn1 = fx₁
    return (state, true)
end

# Main algorithm for A42 method
function update_state(M::A42, F, state::A42State{T,S}, options, l=NullTracks()) where {T,S}
    a::T, b::T, d::T, ee::T = state.xn0, state.xn1, state.d, state.ee
    fa::S, fb::S, fd::S, fe::S = state.fxn0, state.fxn1, state.fd, state.fee

    an, bn = a, b
    μ, λ = 0.5, 0.7
    tole = max(options.xabstol, max(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
    delta = λ * tole

    if isnan(ee)
        c = newton_quadratic(a, b, d, fa, fb, fd, 2)
    else
        c = ipzero(a, b, d, ee, fa, fb, fd, fe)
    end
    fc::S = F(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return _set_state(state, c, fc)
    (isnan(c) || isinf(c)) && return (state, true)

    ab::T, bb::T, db::T, fab::S, fbb::S, fdb::S = bracket(a, b, c, fa, fb, fc)
    eb::T, feb::S = d, fd

    cb::T = take_a42_step(ab, bb, db, eb, fab, fbb, fdb, feb, delta)
    fcb::S = F(cb)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return _set_state(state, c, fc)
    if isnan(c) || isinf(c)
        # tighten up bracket
        @set! state.xn0 = ab
        @set! state.xn1 = bb
        @set! state.d = db
        @set! state.fxn0 = fab
        @set! state.fxn1 = fbb
        @set! state.fd = fdb
        return state, false
    end

    ab, bb, db, fab, fbb, fdb = bracket(ab, bb, cb, fab, fbb, fcb)

    u::T, fu::S = choose_smallest(ab, bb, fab, fbb)
    cb = u - 2 * fu * (bb - ab) / (fbb - fab)
    ch::T = cb
    if abs(cb - u) > 0.5 * (b - a)
        ch = _middle(an, bn)
    end
    fch::S = F(ch)
    incfn(l)

    (iszero(fch) || isnan(fch)) && return _set_state(state, ch, fch)
    if isnan(ch) || isinf(ch)
        # tighten up bracket
        @set! state.xn0 = ab
        @set! state.xn1 = bb
        @set! state.d = db
        @set! state.fxn0 = fab
        @set! state.fxn1 = fbb
        @set! state.fd = fdb
        return state, false
    end

    ah::T, bh::T, dh::T, fah::S, fbh::S, fdh::S = bracket(ab, bb, ch, fab, fbb, fch)

    if bh - ah < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, ee = ah, bh, dh, db
        fa, fb, fd, fe = fah, fbh, fdh, fdb
    else
        m::T = _middle(ah, bh)
        fm::S = F(m)
        incfn(l)
        ee, fe = dh, fdh
        a, b, d, fa, fb, fd = bracket(ah, bh, m, fah, fbh, fm)
    end

    @set! state.xn0 = a
    @set! state.xn1 = b
    @set! state.d = d
    @set! state.ee = ee
    @set! state.fxn0 = fa
    @set! state.fxn1 = fb
    @set! state.fd = fd
    @set! state.fee = fe

    return state, false
end

## ----

"""
    Roots.AlefeldPotraShi()

Follows algorithm in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
EQUATIONS", by Alefeld, Potra, Shi; DOI:
[10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).
Asymptotic efficiency index is ``1.618``. Less efficient, but can run faster than the [`A42`](@ref) method.
Originally by John Travers.
"""
struct AlefeldPotraShi <: AbstractAlefeldPotraShi end

struct AlefeldPotraShiState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    d::T
    fxn1::S
    fxn0::S
    fd::S
end

function init_state(::AlefeldPotraShi, F, x₀, x₁, fx₀, fx₁; c=_middle(x₀, x₁), fc=F(c))
    a, b, fa, fb = x₀, x₁, fx₀, fx₁
    if a > b
        a, b, fa, fb = b, a, fb, fa
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
    return AlefeldPotraShiState(b, a, d, fb, fa, fd)
end
initial_fncalls(::AlefeldPotraShiState) = 3 # worst case assuming fx₀, fx₁,fc must be computed

# ## 3, maybe 4, functions calls per step
function update_state(
    M::AlefeldPotraShi,
    f,
    state::AbstractUnivariateZeroState{T,S},
    options::UnivariateZeroOptions,
    l=NullTracks(),
) where {T,S}
    a::T, b::T, d::T = state.xn0, state.xn1, state.d

    fa::S, fb::S, fd::S = state.fxn0, state.fxn1, state.fd
    μ, λ = 0.5, 0.7
    tole = max(options.xabstol, max(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
    delta = λ * tole

    c::T = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
    fc::S = f(c)
    incfn(l)
    if iszero(fc) || isnan(fc)
        @set! state.xn1 = c
        @set! state.fxn1 = fc
        return (state, true)
    elseif isnan(c) || isinf(c)
        return (state, true)
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    c = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
    fc = f(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return _set_state(state, c, fc)
    if isnan(c) || isinf(c)
        # tighten up bracket
        @set! state.xn0 = a
        @set! state.xn1 = b
        @set! state.d = d
        @set! state.fxn0 = fa
        @set! state.fxn1 = fb
        @set! state.fd = fd

        return (state, false)
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    u::T, fu::S = choose_smallest(a, b, fa, fb)
    c = u - 2 * fu * (b - a) / (fb - fa)
    if abs(c - u) > 0.5 * (b - a)
        c = _middle(a, b)
    end
    fc = f(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return _set_state(state, c, fc)
    if isnan(c) || isinf(c)
        # tighten up bracket
        @set! state.xn0 = a
        @set! state.xn1 = b
        @set! state.d = d
        @set! state.fxn0 = fa
        @set! state.fxn1 = fb
        @set! state.fd = fd

        return (state, false)
    end

    ahat::T, bhat::T, dhat::T, fahat::S, fbhat::S, fdhat::S = bracket(a, b, c, fa, fb, fc)
    if bhat - ahat < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, fa, fb, fd = ahat, bhat, dhat, fahat, fbhat, fdhat
    else
        m::T = _middle(ahat, bhat)
        fm::S = f(m)
        incfn(l)
        a, b, d, fa, fb, fd = bracket(ahat, bhat, m, fahat, fbhat, fm)
    end

    @set! state.xn0 = a
    @set! state.xn1 = b
    @set! state.d = d
    @set! state.fxn0 = fa
    @set! state.fxn1 = fb
    @set! state.fd = fd

    return (state, false)
end

### Brent
"""
    Roots.Brent()

An implementation of
[Brent's](https://en.wikipedia.org/wiki/Brent%27s_method) (or Brent-Dekker) method.
This method uses a choice of inverse quadratic interpolation or a secant
step, falling back on bisection if necessary.

"""
struct Brent <: AbstractBracketing end

struct BrentState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    c::T
    d::T
    fxn1::S
    fxn0::S
    fc::S
    mflag::Bool
end

function log_step(l::Tracks, M::Brent, state)
    a, b = state.xn0, state.xn1
    u, v = a < b ? (a, b) : (b, a)
    push!(l.xs, a)
    push!(l.xs, b) # we store [ai,bi, ai+1, bi+1, ...]
    log_steps(l)
end

# # we store mflag as -1, or +1 in state.mflag
function init_state(::Brent, F, x₀, x₁, fx₀, fx₁)
    u, v, fu, fv = x₀, x₁, fx₀, fx₁
    if abs(fu) > abs(fv)
        u, v, fu, fv = v, u, fv, fu
    end

    BrentState(u, v, v, v, fu, fv, fv, true)
end

function update_state(
    ::Brent,
    f,
    state::BrentState,
    options::UnivariateZeroOptions{T,S},
    l=NullTracks(),
) where {T,S}
    mflag = state.mflag
    a, b, c, d = state.xn0, state.xn1, state.c, state.d
    fa, fb, fc = state.fxn0, state.fxn1, state.fc

    # next step depends on points; inverse quadratic
    s::T = inverse_quadratic_step(a, b, c, fa, fb, fc)
    (isnan(s) || isinf(s)) && (s = secant_step(a, b, fa, fb))

    # guard step
    u, v = (3a + b) / 4, b
    if u > v
        u, v = v, u
    end

    tol = max(options.xabstol, max(abs(b), abs(c), abs(d)) * options.xreltol)
    if !(u < s < v) ||
       (mflag && abs(s - b) >= abs(b - c) / 2) ||
       (!mflag && abs(s - b) >= abs(b - c) / 2) ||
       (mflag && abs(b - c) <= tol) ||
       (!mflag && abs(c - d) <= tol)
        s = _middle(a, b)

        mflag = true
    else
        mflag = false
    end

    fs = f(s)
    incfn(l)

    if iszero(fs)
        @set! state.xn1 = s
        @set! state.fxn1 = fs
        return (state, true)
    elseif isnan(fs) || isinf(fs)
        return (state, true)
    end

    d = c
    c, fc = b, fb

    if sign(fa) * sign(fs) < 0
        b, fb = s, fs
    else
        a, fa = s, fs
    end

    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end

    @set! state.xn0 = a
    @set! state.xn1 = b
    @set! state.c = c
    @set! state.d = d
    @set! state.fxn0 = fa
    @set! state.fxn1 = fb
    @set! state.fc = fc
    @set! state.mflag = mflag

    return state, false
end

## ----------------------------

struct FalsePosition{R} <: AbstractBisection end
"""

    FalsePosition()

Use the [false
position](https://en.wikipedia.org/wiki/False_position_method) method
to find a zero for the function `f` within the bracketing interval
`[a,b]`.

The false position method is a modified bisection method, where the
midpoint between `[aₖ, bₖ]` is chosen to be the intersection point
of the secant line with the ``x`` axis, and not the average between the
two values.

To speed up convergence for concave functions, this algorithm
implements the ``12`` reduction factors of Galdino (*A family of regula
falsi root-finding methods*). These are specified by number, as in
`FalsePosition(2)` or by one of three names `FalsePosition(:pegasus)`,
`FalsePosition(:illinois)`, or `FalsePosition(:anderson_bjork)` (the
default). The default choice has generally better performance than the
others, though there are exceptions.

For some problems, the number of function calls can be greater than
for the `BisectionExact` method, but generally this algorithm will make
fewer function calls.

Examples
```
find_zero(x -> x^5 - x - 1, (-2, 2), FalsePosition())
```
"""
FalsePosition
FalsePosition(x=:anderson_bjork) = FalsePosition{x}()

function default_tolerances(::FalsePosition, ::Type{T}, ::Type{S}) where {T,S}
    default_tolerances(Secant(), T, S)
end

# use fallback for derivative free
function assess_convergence(::FalsePosition, state::UnivariateZeroState, options)
    assess_convergence(Any, state, options)
end

function decide_convergence(
    ::FalsePosition,
    F,
    state::AbstractUnivariateZeroState{T,S},
    options,
    val,
) where {T,S}
    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    isnan(fa) && return b
    isnan(fb) && return a

    if abs(fa) < abs(fb)
        _is_f_approx_0(fa, a, options.abstol, options.reltol, true) && return a
    else
        _is_f_approx_0(fb, b, options.abstol, options.reltol, true) && return b
    end
    val == :not_converged && return NaN * a

    if abs(fa) < abs(fb)
        return a
    else
        return b
    end
end

function update_state(
    method::FalsePosition,
    fs,
    o::AbstractUnivariateZeroState{T,S},
    options::UnivariateZeroOptions,
    l=NullTracks(),
) where {T,S}
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1

    lambda = fb / (fb - fa)
    tau = 1e-10                   # some engineering to avoid short moves; still fails on some
    if !(tau < abs(lambda) < 1 - tau)
        lambda = 1 / 2
    end

    x::T = b - lambda * (b - a)
    fx = fs(x)
    incfn(l)

    if iszero(fx)
        @set! o.xn1 = x
        @set! o.fxn1 = fx
        return (o, true)
    end

    if sign(fx) * sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino_reduction(method, fa, fb, fx)
    end
    b, fb = x, fx

    @set! o.xn0 = a
    @set! o.xn1 = b
    @set! o.fxn0 = fa
    @set! o.fxn1 = fb

    return (o, false)
end

# the 12 reduction factors offered by Galdino
# In RootsTesting.jl, we can see :12 has many more failures.
galdino = Dict{Union{Int,Symbol},Function}(
    :1 => (fa, fb, fx) -> fa * fb / (fb + fx),
    :2 => (fa, fb, fx) -> (fa - fb) / 2,
    :3 => (fa, fb, fx) -> (fa - fx) / (2 + fx / fb),
    :4 => (fa, fb, fx) -> (fa - fx) / (1 + fx / fb)^2,
    :5 => (fa, fb, fx) -> (fa - fx) / (1.5 + fx / fb)^2,
    :6 => (fa, fb, fx) -> (fa - fx) / (2 + fx / fb)^2,
    :7 => (fa, fb, fx) -> (fa + fx) / (2 + fx / fb)^2,
    :8 => (fa, fb, fx) -> fa / 2,
    :9 => (fa, fb, fx) -> fa / (1 + fx / fb)^2,
    :10 => (fa, fb, fx) -> (fa - fx) / 4,
    :11 => (fa, fb, fx) -> fx * fa / (fb + fx),
    :12 => (fa, fb, fx) -> (fa * (1 - fx / fb > 0 ? 1 - fx / fb : 1 / 2)),
)

# give common names
for (nm, i) in [(:pegasus, 1), (:illinois, 8), (:anderson_bjork, 12)]
    galdino[nm] = galdino[i]
end

# from Chris Elrod; https://raw.githubusercontent.com/chriselrod/AsymptoticPosteriors.jl/master/src/false_position.jl
@generated function galdino_reduction(methods::FalsePosition{R}, fa, fb, fx) where {R}
    f = galdino[R]
    quote
        $Expr(:meta, :inline)
        $f(fa, fb, fx)
    end
end

# deprecate this
"""
    find_bracket(f, x0, method=A42(); kwargs...)

For bracketing methods returns an approximate root, the last bracketing interval used, and a flag indicating if an exact zero was found as a named tuple.

With the default tolerances, one  of these should be the case: `exact`  is `true` (indicating termination  of the algorithm due to an exact zero  being identified) or the length of `bracket` is less or equal than `2eps(maximum(abs.(bracket)))`. In the `BisectionExact` case, the 2 could  be replaced by 1, as the bracket, `(a,b)` will satisfy  `nextfloat(a) == b `;  the Alefeld,  Potra, and Shi algorithms don't quite have  that promise.

"""
function find_bracket(
    fs,
    x0,
    method::M=A42();
    kwargs...,
) where {M<:Union{AbstractAlefeldPotraShi,BisectionExact}}
    x = adjust_bracket(x0)
    F = Callable_Function(method, fs) #callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)
    l = NullTracks()

    # check if tolerances are exactly 0
    iszero_tol =
        iszero(options.xabstol) &&
        iszero(options.xreltol) &&
        iszero(options.abstol) &&
        iszero(options.reltol)

    val, stopped = :not_converged, false
    while !stopped
        val, stopped = assess_convergence(method, state, options)
        stopped && break
        state, stopped = update_state(method, F, state, options, l)
    end

    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1
    xstar, fxstar = abs(fa) < abs(fb) ? (a, fa) : (b, fb)
    (xstar=xstar, bracket=(a, b), exact=iszero(fxstar))
end
