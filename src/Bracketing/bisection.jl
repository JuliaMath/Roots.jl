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

initial_fncalls(::Roots.AbstractBisection) = 3

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

function log_step(l::Tracks, M::Bisection, state; init::Bool=false)
    a, b = state.xn0, state.xn1
    push!(l.abₛ, (a,b))
    init && log_iteration(l, 1) # c is computed
    !init && log_iteration(l, 1)
    nothing
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
##################################################
function update_state(M::AbstractBisection, F, o, options, l=NullTracks())
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
