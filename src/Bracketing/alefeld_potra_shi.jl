## --------------------------------------------------

## Two algorithms of Alefeld, Potra, and Shi

## --------------------------------------------------

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

# return c in (a+delta, b-delta)
# adds part of `bracket` from paper with `delta`
function newton_quadratic(a, b, d, fa, fb, fd, k::Int, delta=zero(a))
    a′, b′ = a + 2delta, b - 2delta

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

        if a′ < r < b′
            return r
        end
    end
    # try secant step
    r = secant_step(a, b, fa, fb)

    if a′ < r < b′
        return r
    end

    return _middle(a, b)
end

# f[a, b] divided differences
@inline f_ab(a, b, fa, fb) = (fb - fa) / (b - a)

# f[a,b,d]
@inline function f_abd(a, b, d, fa, fb, fd)
    fab, fbd = f_ab(a, b, fa, fb), f_ab(b, d, fb, fd)
    (fbd - fab) / (d - a)
end

## --------------------------------------------------
##
## Alefeld, Potra, Shi have two algorithms belosw, one is most efficient, but
## slightly slower than other.

abstract type AbstractAlefeldPotraShi <: AbstractBracketingMethod end
initial_fncalls(::AbstractAlefeldPotraShi) = 3 # worst case assuming fx₀, fx₁,fc must be computed

## ----

"""
    Roots.AlefeldPotraShi()

Follows algorithm in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
EQUATIONS", by Alefeld, Potra, Shi; DOI:
[10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).

The order of convergence is `2 + √5`; asymptotically there are 3 function evaluations per step.
Asymptotic efficiency index is ``(2+√5)^(1/3) ≈ 1.618...``. Less efficient, but can run faster than the [`A42`](@ref) method.

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
    isinf(a) && (a = nextfloat(a))
    isinf(b) && (b = prevfloat(b))

    if a > b
        a, b, fa, fb = b, a, fb, fa
    end

    # check if fa*fb ≥ 0
    (iszero(fa) || iszero(fb)) && return AlefeldPotraShiState(b, a, a, fb, fa, fa)
    assert_bracket(fa, fb)

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
    sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

    return AlefeldPotraShiState(b, a, d, fb, fa, fd)
end

# ## 3, maybe 4, functions calls per step
function update_state(
    M::AlefeldPotraShi,
    f,
    state::AlefeldPotraShiState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    a::T, b::T, d::T = state.xn0, state.xn1, state.d

    fa::S, fb::S, fd::S = state.fxn0, state.fxn1, state.fd
    μ, λ = 0.5, 0.7

    tole = max(options.xabstol, min(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
    delta = λ * tole

    c::T = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
    fc::S = f(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
    (isnan(c) || isinf(c)) && return (state, true)

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    c = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
    fc = f(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
    if isnan(c) || isinf(c)
        # tighten up bracket
        state = _set(state, (b, fb), (a, fa))
        @set! state.d = d
        @set! state.fd = fd

        return (state, false)
    end

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    u::T, fu::S = choose_smallest(a, b, fa, fb)
    c = u - 2 * fu * (b - a) / (fb - fa)
    if abs(c - u) > 0.5 * (b - a)
        c = __middle(a, b)
    end
    fc = f(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
    if isnan(c) || isinf(c)
        # tighten up bracket
        state = _set(state, (b, fb), (a, fa))
        @set! state.d = d
        @set! state.fd = fd

        return (state, false)
    end

    ahat::T, bhat::T, dhat::T, fahat::S, fbhat::S, fdhat::S = bracket(a, b, c, fa, fb, fc)
    if bhat - ahat < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, fa, fb, fd = ahat, bhat, dhat, fahat, fbhat, fdhat
    else
        m::T = __middle(ahat, bhat)
        fm::S = f(m)
        incfn(l)
        a, b, d, fa, fb, fd = bracket(ahat, bhat, m, fahat, fbhat, fm)
    end

    state = _set(state, (b, fb), (a, fa))
    @set! state.d = d
    @set! state.fd = fd

    return (state, false)
end

## --------------------------------------------------

# inverse cubic interporation
function ipzero(a, b, c, d, fa, fb, fc, fd)
    Q11 = (c - d) * fc / (fd - fc)
    Q21 = (b - c) * fb / (fc - fb)
    Q31 = (a - b) * fa / (fb - fa)
    D21 = (b - c) * fc / (fc - fb)
    D31 = (a - b) * fb / (fb - fa)
    Q22 = (D21 - Q11) * fb / (fd - fb)
    Q32 = (D31 - Q21) * fa / (fc - fa)
    D32 = (D31 - Q21) * fc / (fc - fa)
    Q33 = (D32 - Q22) * fa / (fd - fa)
    a + (Q31 + Q32 + Q33)
end

# Cubic if possible, if not newton_quadratic until a value is found
function inverse_cubic_interpolation(a, b, d, ee, fa, fb, fd, fe, k, delta=zero(a))
    # if r is NaN or Inf we move on by condition. Faster than checking ahead of time for
    # distinctness
    r = ipzero(a, b, d, ee, fa, fb, fd, fe)
    (a + 2delta < r < b - 2delta) && return r
    r = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
end

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

## initial step, needs to log a,b,d
function log_step(l::Tracks, M::AbstractAlefeldPotraShi, state; init::Bool=false)
    a, b, c = state.xn0, state.xn1, state.d
    init && push!(l.abₛ, extrema((a, b, c)))
    init && log_iteration(l, 1) # take an initial step
    push!(l.abₛ, (a, b))
    !init && log_iteration(l, 1)
    nothing
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
    isinf(a) && (a = nextfloat(a))
    isinf(b) && (b = prevfloat(b))

    if a > b
        a, b, fa, fb = b, a, fb, fa
    end

    # check if fa*fb ≥ 0
    (iszero(fa) || iszero(fb)) && return A42State(b, a, a, a, fb, fa, fa, fa)
    assert_bracket(fa, fb)

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    T = typeof(d)
    ee, fe = T(NaN) / oneunit(T(NaN)) * d, fd # use NaN for initial

    sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

    A42State(b, a, d, ee, fb, fa, fd, fe)
end

# Main algorithm for A42 method
function update_state(M::A42, F, state::A42State{T,S}, options, l=NullTracks()) where {T,S}
    a::T, b::T, d::T, ee::T = state.xn0, state.xn1, state.d, state.ee
    fa::S, fb::S, fd::S, fe::S = state.fxn0, state.fxn1, state.fd, state.fee

    an, bn = a, b
    μ, λ = 0.5, 0.7

    tole = max(options.xabstol, min(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol, not max
    delta = λ * tole

    if isnan(ee)
        c = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
    else
        c = inverse_cubic_interpolation(a, b, d, ee, fa, fb, fd, fe, delta)
    end
    fc::S = F(c)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)

    (isnan(c) || isinf(c)) && return (state, true)

    ab::T, bb::T, db::T, fab::S, fbb::S, fdb::S = bracket(a, b, c, fa, fb, fc)
    eb::T, feb::S = d, fd

    cb::T = inverse_cubic_interpolation(ab, bb, db, eb, fab, fbb, fdb, feb, delta)
    fcb::S = F(cb)
    incfn(l)

    (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
    if isnan(c) || isinf(c)
        # tighten up bracket
        state = _set(state, (bb, fbb), (ab, fab))
        @set! state.d = db
        @set! state.fd = fdb
        return state, false
    end

    ab, bb, db, fab, fbb, fdb = bracket(ab, bb, cb, fab, fbb, fcb)

    u::T, fu::S = choose_smallest(ab, bb, fab, fbb)
    cb = u - 2 * fu * (bb - ab) / (fbb - fab)
    ch::T = cb
    if abs(cb - u) > 0.5 * (b - a)
        ch = __middle(an, bn)
    end
    fch::S = F(ch)
    incfn(l)

    (iszero(fch) || isnan(fch)) && return (_set(state, (ch, fch)), true)
    if isnan(ch) || isinf(ch)
        # tighten up bracket
        state = _set(state, (bb, fbb), (ab, fab))
        @set! state.d = db
        @set! state.fd = fdb
        return state, false
    end

    ah::T, bh::T, dh::T, fah::S, fbh::S, fdh::S = bracket(ab, bb, ch, fab, fbb, fch)

    if bh - ah < μ * (b - a)
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
        a, b, d, ee = ah, bh, dh, db
        fa, fb, fd, fe = fah, fbh, fdh, fdb
    else
        m::T = __middle(ah, bh)
        fm::S = F(m)
        incfn(l)
        ee, fe = dh, fdh
        a, b, d, fa, fb, fd = bracket(ah, bh, m, fah, fbh, fm)
    end

    state = _set(state, (b, fb), (a, fa))
    @set! state.d = d
    @set! state.ee = ee
    @set! state.fd = fd
    @set! state.fee = fe

    return state, false
end
