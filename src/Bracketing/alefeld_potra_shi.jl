## --------------------------------------------------

## Algorithms of Alefeld, Potra, and Shi

## --------------------------------------------------

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
function init_state(::AbstractAlefeldPotraShi, F, x₀, x₁, fx₀, fx₁;
                    #c = secant_step(x₀,x₁,fx₀,fx₁),
                    c = _middle(x₀,x₁)
                    )

    a, b, fa, fb = x₀, x₁, fx₀, fx₁
    iszero(fa) && return AbstractAlefeldPotraShiState(b, a, a, a, fb, fa, fa, fa)
    iszero(fb) && return AbstractAlefeldPotraShiState(b, a, a, a, fb, fa, fa, fa)
#    isinf(a) && (a = nextfloat(a); fa = first(F(a)))
#    isinf(b) && (b = prevfloat(b); fb = first(F(b)))

    if a > b
        a, b, fa, fb = b, a, fb, fa
    end


    c = clamp(c, a, b) # should be unnecessary, but can be with secant step
    fc= first(F(c))

    (iszero(fc) || isnan(fc)) && return AbstractAlefeldPotraShiState(c, a, a, a, fc, fa, fa, fa)

    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
    assert_bracket(fa, fb)

    T = typeof(d)
    ee, fe = T(NaN) / oneunit(T(NaN)) * d, fd # use NaN for initial ee value

    AbstractAlefeldPotraShiState(promote(b, a, d, ee)..., promote(fb, fa, fd, fe)...)
end

# fn calls w/in calculateΔ
# all below use just 1, but others might not (E.g. 4.2 uses K-1
fncalls_per_step(::AbstractAlefeldPotraShi) = 1

function update_state(M::AbstractAlefeldPotraShi,
    F::Callable_Function,
    o::AbstractAlefeldPotraShiState{T,S},
    options,
    l=NullTracks(),
) where {T,S}

    μ, λ = 0.5, 0.7
    atol, rtol = options.xabstol, options.xreltol
    tols = (;λ = λ, atol=atol, rtol=rtol)

    a::T, b::T, d::T, ee::T = o.xn0, o.xn1, o.d, o.ee
    fa::S, fb::S, fd::S, fee::S = o.fxn0, o.fxn1, o.fd, o.fee

    δ₀ = b - a

    # use c to track smaller of |fa|, |fb|
    c, fc = abs(fa) < abs(fb) ? (a,fa) : (b,fb)

    ps = (;a=a, b=b, d=d, ee=ee,
          fa=fa, fb=fb, fd=fd, fee=fee,
          atol=atol, rtol=rtol)

    # (may modify ps) <<----
    Δ, ps = calculateΔ(M, F, c, ps); incfn(l, fncalls_per_step(M))

    a, b, d = ps.a, ps.b, ps.d
    fa, fb, fd = ps.fa, ps.fb, ps.fd

    if iszero(fa) || iszero(fb) || (b - a) <= tolₑ(a,b,fa,fb,atol,rtol)
        @set! o.xn0=a; @set! o.xn1=b; @set! o.fxn0=fa; @set! o.fxn1=fb
        return o, true
   end

    x = c - Δ

    x = avoid_boundaries(a,x,b,fa,fb, tols)
    fx = first(F(x)); incfn(l)

    ā, b̄, d, fā, fb̄, fd = bracket(a,b,x,fa,fb,fx)

    if iszero(fx) || (b̄ - ā) <= tolₑ(ā,b̄,fā,fb̄,atol,rtol)
        @set! o.xn0=ā; @set! o.xn1=b̄; @set! o.fxn0=fā; @set! o.fxn1=fb̄
        return o, true
    end


    u, fu = abs(fā) < abs(fb̄) ? (ā, fā) : (b̄, fb̄)

    # 4.16 double secant step
    fab⁻¹ = (b̄ - ā) / (fb̄ - fā)
    c̄ = u - 2 * fab⁻¹ * fu

    if 2abs(u - c̄) > b̄ - ā
        c̄ = __middle(ā, b̄) #
        #c̄ = ā/2 + b̄/2
    end
    c̄ = avoid_boundaries(ā,c̄,b̄,fā,fb̄,  tols)
    fc̄ = first(F(c̄)); incfn(l)
    (iszero(fc̄) ||isnan(fc̄)) && return (_set(o, (c̄, fc̄)), true)

    â,b̂,d̂,fâ,fb̂,fd̂ = bracket(ā,b̄,c̄,fā,fb̄,fc̄)
    ##
    if (b̂ - â) < μ*δ₀
        ee, fee = d, fd
        a, b, d, fa, fb, fd = â, b̂, d̂, fâ, fb̂, fd̂
    else
        m =  __middle(ā, b̄)
        #m = ā/2 + b̄/2
        m = avoid_boundaries(â, m, b̂, fâ, fb̂, tols)
        fm = first(F(m)); incfn(l)
        (iszero(fm)||isnan(fm)) && return (_set(o, (m, fm)), true)

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

struct A2425{K} <: AbstractAlefeldPotraShi end
function calculateΔ(::A2425{K}, F::Callable_Function, c₀::T, ps) where {K, T}
    a,b,d,ee = ps.a, ps.b, ps.d, ps.ee
    fa,fb,fd,fee = ps.fa, ps.fb, ps.fd, ps.fee
    tols = (λ=0.7, atol=ps.atol, rtol=ps.rtol)

    c = a
    for k ∈ 1:K
        c = newton_quadratic(a, b, d, fa, fb, fd, k+1)

        k == K && break
        c = avoid_boundaries(a,c,b,fa,fb,  tols)
        fc = first(F(c))
        a,b,d,fa,fb,fd = bracket(a,b,c,fa,fb,fc)

        iszero(fc) && break
        if (isnan(fc) || isnan(c) || isinf(c))
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

    c₀ - c,  ps


end


"""
    Roots.AlefeldPotraShi()

Follows algorithm in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
EQUATIONS", by Alefeld, Potra, Shi; DOI:
[10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).

The order of convergence is `2 + √5`; asymptotically there are 3 function evaluations per step.
Asymptotic efficiency index is ``(2+√5)^{1/3} ≈ 1.618...``. Less efficient, but can run faster than the [`A42`](@ref) method.

Originally by John Travers.
"""
const AlefeldPotraShi = A2425{2}


# Algorith 5.7 is parameterized by K
# 4.1 -> K=1; 4.2 -> K=2
struct A57{K} <: AbstractAlefeldPotraShi end
fncalls_per_step(::A57{K}) where {K} = K - 1
function calculateΔ(::A57{K}, F::Callable_Function, c₀::T, ps) where {K, T}
    a,b,d,ee = ps.a, ps.b, ps.d, ps.ee
    fa,fb,fd,fee = ps.fa, ps.fb, ps.fd, ps.fee
    tols = (λ=0.7, atol=ps.atol, rtol=ps.rtol)
    c, fc = a, fa

    for k ∈ 1:K
        if isnan(ee) || iszero(_pairwise_prod(fa,fb,fd,fee))
            c = newton_quadratic(a, b, d, fa, fb, fd, k+1)
        else
            c = ipzero(a, b, d, ee, fa, fb, fd, fee)
            if (c <= a || b <= c)
                c = newton_quadratic(a, b, d, fa, fb, fd, k+1)
            end
        end

        k == K && break

        ee, fee = d, fd
        c = avoid_boundaries(a, c, b, fa, fb,  tols)
        fc = first(F(c))
        a,b,d,fa,fb,fd = bracket(a,b,c,fa,fb,fc)

        iszero(fc) && break # fa or fb is 0
        if (isnan(fc) || isnan(c) || isinf(c))
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

    c₀ - c,  ps


end

# """
#     Roots.A42()

# Bracketing method which finds the root of a continuous function within
# a provided interval `[a, b]`, without requiring derivatives. It is based
# on algorithm 4.2 described in: G. E. Alefeld, F. A. Potra, and
# Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
# Trans. Math. Softw. 21, 327–344 (1995), DOI: [10.1145/210089.210111](https://doi.org/10.1145/210089.210111).
# The asymptotic efficiency index, ``q^{1/k}``, is ``(2 + 7^{1/2})^{1/3} = 1.6686...``.


# Originally by John Travers.

# """
const A42 = A57{2}


## --- utilities

## Brent-style tole from paper
function tolₑ(a,b,fa,fb, atol,rtol)
    u = abs(fa) < abs(fb) ? abs(a) : abs(b)
    return 2*u*rtol + atol
end

## adjustment before calling bracket
function avoid_boundaries(a,c,b, fa, fb, tols)

    δ = tols.λ * tolₑ(a, b, fa, fb, tols.atol, tols.rtol)

    if (b - a) ≤ 4δ
        c = a/2 + b/2
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


# different from above; no δ adjustment
function newton_quadratic(a, b, d, fa, fb, fd, k::Int)

    A = f_abd(a, b, d, fa, fb, fd)
    B = f_ab(a, b, fa, fb)

    (iszero(A) || isinf(A) || isnan(A)) && return a - fa / B


    r = sign(A)*sign(fa) > 0 ? a : b

    for i in 1:k
        P = fa + B * (r - a) + A * (r - a) * (r - b)
        P′ = (B + A * (2r - a - b))
        r -= P / P′
    end

    return r
end

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

function _pairwise_prod(as...)
    t = one(first(as))
    n = length(as)
    for i ∈ 1:n-1
        for j∈ (i+1):n
            t *= (as[i]-as[j])
        end
    end
    t
end

## ------- delete me ------

## --- older versions. Delete if no issues with newer
## try to make a general framework for alefeld potra shi
## *If* this is as efficient as A42 and AlefeldPotraShi above, we
## can remove them, for now, this isn't
## can opt in by creating calculateΔ function for a method
## points a,b,d,e stored

# return c in (a+delta, b-delta)
# adds part of `bracket` from paper with `delta`
# function newton_quadratic(a, b, d, fa, fb, fd, k::Int, delta=zero(a))
#     a′, b′ = a + 2delta, b - 2delta

#     A = f_abd(a, b, d, fa, fb, fd)

#     r = isbracket(A, fa) ? b : a

#     # use quadratic step; if that fails, use secant step; if that fails, bisection
#     if !(isnan(A) || isinf(A)) || !iszero(A)
#         B = f_ab(a, b, fa, fb)

#         for i in 1:k
#             Pr = fa + B * (r - a) + A * (r - a) * (r - b)
#             Prp = (B + A * (2r - a - b))
#             r -= Pr / Prp
#         end

#         if a′ < r < b′
#             return r
#         end
#     end
#     # try secant step
#     r = secant_step(a, b, fa, fb)

#     if a′ < r < b′
#         return r
#     end

#     return _middle(a, b)
# end

# # inverse cubic interporation

# # Cubic if possible, if not newton_quadratic until a value is found
# function inverse_cubic_interpolation(a, b, d, ee, fa, fb, fd, fe, k, delta=zero(a))
#     # if r is NaN or Inf we move on by condition. Faster than checking ahead of time for
#     # distinctness
#     r = ipzero(a, b, d, ee, fa, fb, fd, fe)
#     (a + 2delta < r < b - 2delta) && return r
#     r = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
# end



## --------------------------------------------------
##
## Alefeld, Potra, Shi have two algorithms belosw, one is most efficient, but
## slightly slower than other.


# ## ----

# """
#     Roots.AlefeldPotraShi()

# Follows algorithm in "ON ENCLOSING SIMPLE ROOTS OF NONLINEAR
# EQUATIONS", by Alefeld, Potra, Shi; DOI:
# [10.1090/S0025-5718-1993-1192965-2](https://doi.org/10.1090/S0025-5718-1993-1192965-2).

# The order of convergence is `2 + √5`; asymptotically there are 3 function evaluations per step.
# Asymptotic efficiency index is ``(2+√5)^{1/3} ≈ 1.618...``. Less efficient, but can run faster than the [`A42`](@ref) method.

# Originally by John Travers.
# """
# struct AlefeldPotraShi <: AbstractAlefeldPotraShi end

# struct AlefeldPotraShiState{T,S} <: AbstractUnivariateZeroState{T,S}
#     xn1::T
#     xn0::T
#     d::T
#     fxn1::S
#     fxn0::S
#     fd::S
# end

# function init_state(::AlefeldPotraShi, F, x₀, x₁, fx₀, fx₁; c=_middle(x₀, x₁), fc=F(c))
#     a, b, fa, fb = x₀, x₁, fx₀, fx₁
#     isinf(a) && (a = nextfloat(a))
#     isinf(b) && (b = prevfloat(b))

#     if a > b
#         a, b, fa, fb = b, a, fb, fa
#     end

#     # check if fa*fb ≥ 0
#     (iszero(fa) || iszero(fb)) && return AlefeldPotraShiState(b, a, a, fb, fa, fa)
#     assert_bracket(fa, fb)

#     a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
#     sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

#     return AlefeldPotraShiState(promote(b, a, d)..., promote(fb, fa, fd)...)
# end

# # ## 3, maybe 4, functions calls per step
# function update_state(
#     M::AlefeldPotraShi,
#     f,
#     state::AlefeldPotraShiState{T,S},
#     options,
#     l=NullTracks(),
# ) where {T,S}
#     a::T, b::T, d::T = state.xn0, state.xn1, state.d

#     fa::S, fb::S, fd::S = state.fxn0, state.fxn1, state.fd
#     μ, λ = 0.5, 0.7

#     tole = max(options.xabstol, min(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol
#     delta = λ * tole

#     c::T = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
#     fc::S = f(c)
#     incfn(l)

#     (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
#     (isnan(c) || isinf(c)) && return (state, true)

#     a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

#     c = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
#     fc = f(c)
#     incfn(l)

#     (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
#     if isnan(c) || isinf(c)
#         # tighten up bracket
#         state = _set(state, (b, fb), (a, fa))
#         @set! state.d = d
#         @set! state.fd = fd

#         return (state, false)
#     end

#     a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

#     u::T, fu::S = choose_smallest(a, b, fa, fb)
#     c = u - 2 * fu * (b - a) / (fb - fa)
#     if abs(c - u) > 0.5 * (b - a)
#         c = __middle(a, b)
#     end
#     fc = f(c)
#     incfn(l)

#     (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
#     if isnan(c) || isinf(c)
#         # tighten up bracket
#         state = _set(state, (b, fb), (a, fa))
#         @set! state.d = d
#         @set! state.fd = fd

#         return (state, false)
#     end

#     ahat::T, bhat::T, dhat::T, fahat::S, fbhat::S, fdhat::S = bracket(a, b, c, fa, fb, fc)
#     if bhat - ahat < μ * (b - a)
#         #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
#         a, b, d, fa, fb, fd = ahat, bhat, dhat, fahat, fbhat, fdhat
#     else
#         m::T = __middle(ahat, bhat)
#         fm::S = f(m)
#         incfn(l)
#         a, b, d, fa, fb, fd = bracket(ahat, bhat, m, fahat, fbhat, fm)
#     end

#     state = _set(state, (b, fb), (a, fa))
#     @set! state.d = d
#     @set! state.fd = fd

#     return (state, false)
# end

## --------------------------------------------------


# """
#     Roots.A42()

# Bracketing method which finds the root of a continuous function within
# a provided interval `[a, b]`, without requiring derivatives. It is based
# on algorithm 4.2 described in: G. E. Alefeld, F. A. Potra, and
# Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
# Trans. Math. Softw. 21, 327–344 (1995), DOI: [10.1145/210089.210111](https://doi.org/10.1145/210089.210111).
# The asymptotic efficiency index, ``q^{1/k}``, is ``(2 + 7^{1/2})^{1/3} = 1.6686...``.


# Originally by John Travers.

# """
# struct A42 <: AbstractAlefeldPotraShi end

# ## initial step, needs to log a,b,d
# function log_step(l::Tracks, M::AbstractAlefeldPotraShi, state; init::Bool=false)
#     a, b, c = state.xn0, state.xn1, state.d
#     init && push!(l.abₛ, extrema((a, b, c)))
#     init && log_iteration(l, 1) # take an initial step
#     push!(l.abₛ, (a, b))
#     !init && log_iteration(l, 1)
#     nothing
# end

# struct A42State{T,S} <: AbstractUnivariateZeroState{T,S}
#     xn1::T
#     xn0::T
#     d::T
#     ee::T
#     fxn1::S
#     fxn0::S
#     fd::S
#     fee::S
# end

# function init_state(::A42, F, x₀, x₁, fx₀, fx₁; c=_middle(x₀, x₁), fc=F(c))
#     a, b, fa, fb = x₀, x₁, fx₀, fx₁
#     isinf(a) && (a = nextfloat(a))
#     isinf(b) && (b = prevfloat(b))

#     if a > b
#         a, b, fa, fb = b, a, fb, fa
#     end

#     # check if fa*fb ≥ 0
#     (iszero(fa) || iszero(fb)) && return A42State(b, a, a, a, fb, fa, fa, fa)
#     assert_bracket(fa, fb)

#     a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

#     T = typeof(d)
#     ee, fe = T(NaN) / oneunit(T(NaN)) * d, fd # use NaN for initial

#     sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

#     A42State(b, a, d, ee, fb, fa, fd, fe)
# end

# # Main algorithm for A42 method
# function update_state(M::A42, F, state::A42State{T,S}, options, l=NullTracks()) where {T,S}
#     a::T, b::T, d::T, ee::T = state.xn0, state.xn1, state.d, state.ee
#     fa::S, fb::S, fd::S, fe::S = state.fxn0, state.fxn1, state.fd, state.fee

#     an, bn = a, b
#     μ, λ = 0.5, 0.7

#     tole = max(options.xabstol, min(abs(a), abs(b)) * options.xreltol) # paper uses 2|u|*rtol + atol, not max
#     delta = λ * tole

#     if isnan(ee)
#         c = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
#     else
#         c = inverse_cubic_interpolation(a, b, d, ee, fa, fb, fd, fe, delta)
#     end
#     fc::S = F(c)
#     incfn(l)

#     (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)

#     (isnan(c) || isinf(c)) && return (state, true)

#     ab::T, bb::T, db::T, fab::S, fbb::S, fdb::S = bracket(a, b, c, fa, fb, fc)
#     eb::T, feb::S = d, fd

#     cb::T = inverse_cubic_interpolation(ab, bb, db, eb, fab, fbb, fdb, feb, delta)
#     fcb::S = F(cb)
#     incfn(l)

#     (iszero(fc) || isnan(fc)) && return (_set(state, (c, fc)), true)
#     if isnan(c) || isinf(c)
#         # tighten up bracket
#         state = _set(state, (bb, fbb), (ab, fab))
#         @set! state.d = db
#         @set! state.fd = fdb
#         return state, false
#     end

#     ab, bb, db, fab, fbb, fdb = bracket(ab, bb, cb, fab, fbb, fcb)

#     u::T, fu::S = choose_smallest(ab, bb, fab, fbb)
#     cb = u - 2 * fu * (bb - ab) / (fbb - fab)
#     ch::T = cb
#     if abs(cb - u) > 0.5 * (b - a)
#         ch = __middle(an, bn)
#     end
#     fch::S = F(ch)
#     incfn(l)

#     (iszero(fch) || isnan(fch)) && return (_set(state, (ch, fch)), true)
#     if isnan(ch) || isinf(ch)
#         # tighten up bracket
#         state = _set(state, (bb, fbb), (ab, fab))
#         @set! state.d = db
#         @set! state.fd = fdb
#         return state, false
#     end

#     ah::T, bh::T, dh::T, fah::S, fbh::S, fdh::S = bracket(ab, bb, ch, fab, fbb, fch)

#     if bh - ah < μ * (b - a)
#         #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat # typo in paper
#         a, b, d, ee = ah, bh, dh, db
#         fa, fb, fd, fe = fah, fbh, fdh, fdb
#     else
#         m::T = __middle(ah, bh)
#         fm::S = F(m)
#         incfn(l)
#         ee, fe = dh, fdh
#         a, b, d, fa, fb, fd = bracket(ah, bh, m, fah, fbh, fm)
#     end

#     state = _set(state, (b, fb), (a, fa))
#     @set! state.d = d
#     @set! state.ee = ee
#     @set! state.fd = fd
#     @set! state.fee = fe

#     return state, false
# end

## ----

# function default_tolerances(::AbstractAlefeldPotraShi, ::Type{T}, ::Type{S}) where {T,S}
#     xatol = zero(real(T)) * oneunit(real(T))
#     xrtol = eps(real(T))  # unitless
#     atol = 0 * oneunit(real(S))
#     rtol = 0 * one(real(S))
#     maxevals = 60
#     strict = true
#     (xatol, xrtol, atol, rtol, maxevals, strict)
# end





## ----

function newton_quadratic(a, b, d, fa, fb, fd, k::Int, delta)
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


# Cubic if possible, if not newton_quadratic until a value is found
function inverse_cubic_interpolation(a, b, d, ee, fa, fb, fd, fe, k, delta=zero(a))
    # if r is NaN or Inf we move on by condition. Faster than checking ahead of time for
    # distinctness
    r = ipzero(a, b, d, ee, fa, fb, fd, fe)

    (a + 2delta < r < b - 2delta) && return r

    r = newton_quadratic(a, b, d, fa, fb, fd, 3, delta)
end
struct A42O <: AbstractAlefeldPotraShi end


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

function init_state(::A42O, F, x₀, x₁, fx₀, fx₁; c=_middle(x₀, x₁), fc=F(c))
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
function update_state(M::A42O, F, state::A42State{T,S}, options, l=NullTracks()) where {T,S}
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
