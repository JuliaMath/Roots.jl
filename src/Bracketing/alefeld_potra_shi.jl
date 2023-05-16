## --------------------------------------------------

## Algorithms of Alefeld, Potra, and Shi
## Main paper is
## [Algorithm 748: enclosing zeros of continuous functions](https://na.math.kit.edu/alefeld/download/1995_Algorithm_748_Enclosing_Zeros_of_Continuous_Functions.pdf) with supporting material at
## [Algorithm 748: enclosing zeros of continuous functions](https://doi.org/10.1145/210089.210111)

## --------------------------------------------------

"""
    AbstractAlefeldPotraShi

An abstract type for Alefeld-Potra-Shi type bracketing problems, as discussed in  [*Algorithm 748: enclosing zeros of continuous functions*](https://na.math.kit.edu/alefeld/download/1995_Algorithm_748_Enclosing_Zeros_of_Continuous_Functions.pdf). These consist of an interpolation step, such as quadratic interpolation or inverse cubic interpolation along with a double-secant step. For a smooth function and finite bracketing interval, these methods should always converge.

The `update_step` method calls a `calculateΔ` method that can be customized to turn an algorithm based on interpolation into a bracketed algorithm. See [`Roots.BracketedHalley`](@ref) for an example.

This implementation deviates slightly from the printed algorithm, as it may use an initial call to `_middle` rather than a secant step, depending on the signs of ``a`` and ``b``.
"""
abstract type AbstractAlefeldPotraShi <: AbstractBracketingMethod end

initial_fncalls(::AbstractAlefeldPotraShi) = 3 # worst case assuming fx₀, fx₁,fc must be computed

## initial step, needs to log a,b,d
function log_step(l::Tracks, M::AbstractAlefeldPotraShi, state; init::Bool=false)
    a, b, c = state.xn0, state.xn1, state.d
    init && push!(l.abₛ, extrema((a, b, c)))
    init && log_iteration(l, 1) # take an initial step
    push!(l.abₛ, (a, b))
    !init && log_iteration(l, 1)
    nothing
end

struct AbstractAlefeldPotraShiState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    d::T
    ee::T
    fxn1::S
    fxn0::S
    fd::S
    fee::S
end

# basic init state is like bracketing state
# keep a < b
# set d, ee to a
function init_state(::AbstractAlefeldPotraShi, F, x₀, x₁, fx₀, fx₁; c=nothing, fc=nothing)
    a, b, fa, fb = x₀, x₁, fx₀, fx₁
    iszero(fa) && return AbstractAlefeldPotraShiState(b, a, a, a, fb, fa, fa, fa)
    iszero(fb) && return AbstractAlefeldPotraShiState(b, a, a, a, fb, fa, fa, fa)

    if a > b
        a, b, fa, fb = b, a, fb, fa
    end

    if c === nothing # need c, fc to be defined if one is
        c = a < zero(a) < b ? _middle(a, b) : secant_step(a, b, fa, fb)
        fc = first(F(c))
    end

    (iszero(fc) || !isfinite(fc)) &&
        return AbstractAlefeldPotraShiState(c, a, a, a, fc, fa, fa, fa)

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
    assert_bracket(fa, fb)

    T = typeof(d)
    ee, fe = T(NaN) / oneunit(T(NaN)) * d, fd # use NaN for initial ee value

    AbstractAlefeldPotraShiState(promote(b, a, d, ee)..., promote(fb, fa, fd, fe)...)
end

# fn calls w/in calculateΔ
# 1 is default, but this should be adjusted for different methods
fncalls_per_step(::AbstractAlefeldPotraShi) = 1

function update_state(
    M::AbstractAlefeldPotraShi,
    F::Callable_Function,
    o::AbstractAlefeldPotraShiState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    μ, λ = 0.5, 0.7
    atol, rtol = options.xabstol, options.xreltol
    tols = (; λ=λ, atol=atol, rtol=rtol)

    a::T, b::T, d::T, ee::T = o.xn0, o.xn1, o.d, o.ee
    fa::S, fb::S, fd::S, fee::S = o.fxn0, o.fxn1, o.fd, o.fee

    δ₀ = b - a

    # use c to track smaller of |fa|, |fb|
    c, fc = abs(fa) < abs(fb) ? (a, fa) : (b, fb)

    ps = (; a=a, b=b, d=d, ee=ee, fa=fa, fb=fb, fd=fd, fee=fee, atol=atol, rtol=rtol)

    # (may modify ps) <<----
    Δ, ps = calculateΔ(M, F, c, ps)
    incfn(l, fncalls_per_step(M))

    a, b, d = ps.a, ps.b, ps.d
    fa, fb, fd = ps.fa, ps.fb, ps.fd

    if iszero(fa) || iszero(fb) || (b - a) <= tolₑ(a, b, fa, fb, atol, rtol)
        @set! o.xn0 = a
        @set! o.xn1 = b
        @set! o.fxn0 = fa
        @set! o.fxn1 = fb
        return o, true
    end

    x = c - Δ

    x = avoid_boundaries(a, x, b, fa, fb, tols)
    fx = first(F(x))
    incfn(l)

    ā, b̄, d, fā, fb̄, fd = bracket(a, b, x, fa, fb, fx)

    if iszero(fx) || (b̄ - ā) <= tolₑ(ā, b̄, fā, fb̄, atol, rtol)
        @set! o.xn0 = ā
        @set! o.xn1 = b̄
        @set! o.fxn0 = fā
        @set! o.fxn1 = fb̄
        return o, true
    end

    u, fu = abs(fā) < abs(fb̄) ? (ā, fā) : (b̄, fb̄)

    # 4.16 double secant step
    fab⁻¹ = (b̄ - ā) / (fb̄ - fā)
    c̄ = u - 2 * fab⁻¹ * fu                  # 4.16

    if 2abs(u - c̄) > b̄ - ā                  # 4.17
        c̄ = __middle(ā, b̄)
    end

    c̄ = avoid_boundaries(ā, c̄, b̄, fā, fb̄, tols)
    fc̄ = first(F(c̄))
    incfn(l)
    (iszero(fc̄) || !isfinite(fc̄)) && return (_set(o, (c̄, fc̄)), true)

    â, b̂, d̂, fâ, fb̂, fd̂ = bracket(ā, b̄, c̄, fā, fb̄, fc̄) # 4.18

    if (b̂ - â) < μ * δ₀                        # 4.19
        ee, fee = d, fd
        a, b, d, fa, fb, fd = â, b̂, d̂, fâ, fb̂, fd̂
    else
        m = __middle(ā, b̄)
        m = avoid_boundaries(â, m, b̂, fâ, fb̂, tols)

        fm = first(F(m))
        incfn(l)
        (iszero(fm) || !isfinite(fm)) && return (_set(o, (m, fm)), true)

        ee, fee = d̂, fd̂
        a, b, d, fa, fb, fd = bracket(â, b̂, m, fâ, fb̂, fm)
    end

    @set! o.xn0 = a
    @set! o.xn1 = b
    @set! o.d = d
    @set! o.ee = ee
    @set! o.fxn0 = fa
    @set! o.fxn1 = fb
    @set! o.fd = fd
    @set! o.fee = fee

    return o, false
end

## --- Methods
# algorithms 2.4 and 2.5 can be implemented this way:
struct A2425{K} <: AbstractAlefeldPotraShi end
function calculateΔ(::A2425{K}, F::Callable_Function, c₀::T, ps) where {K,T}
    a, b, d, ee = ps.a, ps.b, ps.d, ps.ee
    fa, fb, fd, fee = ps.fa, ps.fb, ps.fd, ps.fee
    tols = (λ=0.7, atol=ps.atol, rtol=ps.rtol)

    c = a
    for k in 1:K
        c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)

        k == K && break
        c = avoid_boundaries(a, c, b, fa, fb, tols)
        fc = first(F(c))
        a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

        iszero(fc) && break
        if (isnan(fc) || !isfinite(c))
            c = c₀
            break
        end
    end

    @set! ps.a = a
    @set! ps.fa = fa
    @set! ps.b = b
    @set! ps.fb = fb
    @set! ps.d = d
    @set! ps.fd = fd

    c₀ - c, ps
end

"""
    Roots.AlefeldPotraShi()

Follows algorithm 4.1 in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
EQUATIONS", by Alefeld, Potra, Shi; DOI:
[10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).

The order of convergence is `2 + √5`; asymptotically there are 3 function evaluations per step.
Asymptotic efficiency index is ``(2+√5)^{1/3} ≈ 1.618...``. Less efficient, but can run faster than the related [`A42`](@ref) method.

Originally by John Travers.
"""
const AlefeldPotraShi = A2425{2}

# Algorith 5.7 is parameterized by K
# 4.1 -> K=1; 4.2 -> K=2
struct A57{K} <: AbstractAlefeldPotraShi end
fncalls_per_step(::A57{K}) where {K} = K - 1
function calculateΔ(::A57{K}, F::Callable_Function, c₀::T, ps) where {K,T}
    a, b, d, ee = ps.a, ps.b, ps.d, ps.ee
    fa, fb, fd, fee = ps.fa, ps.fb, ps.fd, ps.fee
    tols = (λ=0.7, atol=ps.atol, rtol=ps.rtol)
    c, fc = a, fa

    for k in 1:K
        if isnan(ee) || iszero(_pairwise_prod(fa, fb, fd, fee))
            c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
        else
            c = ipzero(a, b, d, ee, fa, fb, fd, fee)
            if (c <= a || b <= c)
                c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
            end
        end

        k == K && break

        ee, fee = d, fd
        c = avoid_boundaries(a, c, b, fa, fb, tols)
        fc = first(F(c))
        a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

        iszero(fc) && break # fa or fb is 0
        if (!isfinite(fc) || !isfinite(c))
            c = c₀
            break
        end
    end
    @set! ps.a = a
    @set! ps.fa = fa
    @set! ps.b = b
    @set! ps.fb = fb
    @set! ps.d = d
    @set! ps.fd = fd
    @set! ps.ee = ee
    @set! ps.fee = fee

    c₀ - c, ps
end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided bracketing interval `[a, b]`, without requiring derivatives. It is based
on algorithm 4.2 described in: G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995), DOI: [10.1145/210089.210111](https://doi.org/10.1145/210089.210111).
The asymptotic efficiency index, ``q^{1/k}``, is ``(2 + 7^{1/2})^{1/3} = 1.6686...``.

Originally by John Travers.

!!! note

    The paper refenced above shows that for a continuously differentiable ``f`` over ``[a,b]`` with a simple root  the algorithm terminates at a zero or asymptotically the steps are of the inverse cubic type (Lemma 5.1). This is proved under an assumption that ``f`` is four-times continuously differentiable.

"""
const A42 = A57{2}

## --- utilities

## Brent-style tole from paper
function tolₑ(a, b, fa, fb, atol, rtol)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    return 2 * u * rtol + atol
end

## adjustment before calling bracket
function avoid_boundaries(a, c, b, fa, fb, tols)
    δ = tols.λ * tolₑ(a, b, fa, fb, tols.atol, tols.rtol)

    if (b - a) ≤ 4δ
        c = a / 2 + b / 2
    elseif c ≤ a + 2δ
        c = a + 2δ
    elseif c ≥ b - 2δ
        c = b - 2δ
    end
    c
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

# f[a, b] divided differences
@inline f_ab(a, b, fa, fb) = (fb - fa) / (b - a)

# f[a,b,d]
@inline function f_abd(a, b, d, fa, fb, fd)
    fab, fbd = f_ab(a, b, fa, fb), f_ab(b, d, fb, fd)
    (fbd - fab) / (d - a)
end

# iterative quadratic solution to P(x) = 0 where P=f(a) + f[a,b]*(x-a) + f[a,b,d]*(x-a)*(x-b)
function newton_quadratic(a, b, d, fa, fb, fd, k::Int)
    A = f_abd(a, b, d, fa, fb, fd)
    B = f_ab(a, b, fa, fb)

    (iszero(A) || !isfinite(A)) && return a - fa / B

    r = sign(A) * sign(fa) > 0 ? a : b

    for i in 1:k
        P = fa + B * (r - a) + A * (r - a) * (r - b)
        P′ = (B + A * (2r - a - b))
        r -= P / P′
    end

    return r
end

# zero of inverse interpolation polynomial through (a,fa), (b,fb), (c,fc), (d, fd)
# may not lie in [a,b], though asymptotically will under smoothness assumptions
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

# check if fa,fb,fc,fd are distinct
function _pairwise_prod(as...)
    t = one(first(as))
    n = length(as)
    for i in 1:(n - 1)
        for j in (i + 1):n
            t *= (as[i] - as[j])
        end
    end
    t
end
