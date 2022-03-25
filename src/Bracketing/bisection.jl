"""

    Bisection()

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
length is less than or equal to the tolerance `max(xtol, min(abs(a),
abs(b)) * xrtol)` or the function value is less than `max(tol,
min(abs(a), abs(b)) * rtol)`. The latter is used only if the default
tolerances are adjusted.


"""
struct Bisection <: AbstractBisectionMethod end


initial_fncalls(::AbstractBisectionMethod) = 3 # middle

# Bisection using __middle should have a,b on same side of 0.0 (though
# possibly, may be -0.0, 1.0 so not guaranteed to be of same sign)
function init_state(::AbstractBisectionMethod, F, x₀, x₁, fx₀, fx₁; m=_middle(x₀, x₁), fm=F(m))

    if x₀ > x₁
        x₀, x₁, fx₀, fx₁ = x₁, x₀, fx₁, fx₀
    end

    # handle interval if fa*fb ≥ 0 (explicit, but also not needed)
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)
    assert_bracket(fx₀, fx₁)

    if sign(fm) * fx₀ < 0 * oneunit(fx₀)
        a, b, fa, fb = x₀, m, fx₀, fm
    else
        a, b, fa, fb = m, x₁, fm, fx₁
    end

    # handles case where a=-0.0, b=1.0 without error
    sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

    UnivariateZeroState(b, a, fb, fa)
end

const FloatNN = Union{Float64,Float32,Float16}

# for Bisection, the defaults are zero tolerances and strict=true
"""
    default_tolerances(M::AbstractBisectionMethod, [T], [S])

For `Bisection` when the `x` values are of type `Float64`, `Float32`,
or `Float16`, the default tolerances are zero and there is no limit on
the number of iterations. In this case, the
algorithm is guaranteed to converge to an exact zero, or a point where
the function changes sign at one of the answer's adjacent floating
point values.

For other types, default non-zero tolerances for `xatol` and `xrtol` are given.

"""
function default_tolerances(::AbstractBisectionMethod, ::Type{T}, ::Type{S′}) where {T<:FloatNN,S′}
    S = real(float(S′))
    xatol = 0 * oneunit(S)
    xrtol = 0 * one(T)
    atol = 0 * oneunit(S)
    rtol = 0 * one(S)
    maxevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, strict)
end

# not float uses some non-zero tolerances for `x`
function default_tolerances(::AbstractBisectionMethod, ::Type{T′}, ::Type{S′}) where {T′,S′}
    T,S = real(float(T′)), real(float(S′))
    xatol = eps(T)^3 * oneunit(T)
    xrtol = eps(T) * one(T) # unitless
    atol = 0 * oneunit(S)
    rtol = 0 * one(S)
    maxevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, strict)
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

## find middle assuming a,b same sign, finite
## Alternative "mean" definition that operates on the binary representation
## of a float. Using this definition, bisection will never take more than
## 64 steps (over Float64)
__middle(x::Float64, y::Float64) = __middle(Float64, UInt64, x, y)
__middle(x::Float32, y::Float32) = __middle(Float32, UInt32, x, y)
__middle(x::Float16, y::Float16) = __middle(Float16, UInt16, x, y)
## fallback for non FloatNN number types
__middle(x::Number, y::Number) = x / 2 + y / 2

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
# function assess_convergence(::Bisection, state::AbstractUnivariateZeroState, options)

#     a, b = state.xn0, state.xn1
#     fa, fb = state.fxn0, state.fxn1

#     m = min(abs(a), abs(b))
#     δₓ = max(options.xabstol, m * options.xreltol)
#     if δₓ/oneunit(δₓ) == 0
#         (iszero(fa) || isnan(fa) || iszero(fb) || isnan(fb)) && return (:f_converged, true)
#         nextfloat(a) == b && return (:x_converged, true)
#     else
#         abs(b-a) ≤ δₓ && return (:x_converged, true)
#     end

#     δ  = max(options.abstol, (m / oneunit(m)) * (options.reltol * oneunit(fa)))

#     if δ/oneunit(δ) != 0
#         min(abs(fa), abs(fb)) ≤ δ && return (:f_converged, true)
#     end


#     return :not_converged, false
# end

## --------------------------------------------------

function update_state(M::AbstractBisectionMethod, F, o, options, l=NullTracks())
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1

    c = __middle(a, b)
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
