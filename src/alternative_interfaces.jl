## Some **legacy** alternative interfaces.

## several named interfaces to methods
## newton, halley, quadratic_inverse, superhalley, chebyshevlike
"""
    Roots.newton(f, fp, x0; kwargs...)

Implementation of Newton's method: `xᵢ₊₁ =  xᵢ - f(xᵢ)/f'(xᵢ)`.

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`.

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero` using the `Roots.Newton()` method.

See also `Roots.newton((f,fp), x0)` and `Roots.newton(fΔf, x0)` for simpler implementations.

"""
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)

## --------------------------------------------------

"""
    Roots.halley(f, fp, fpp, x0; kwargs...)

Implementation of Halley's method (cf `?Roots.Halley()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.Halley()` method.


"""
halley(f, fp, fpp, x0; kwargs...) = find_zero((f, fp, fpp), x0, Halley(); kwargs...)

"""
    Roots.quadratic_inverse(f, fp, fpp, x0; kwargs...)

Implementation of the quadratic inverse method (cf `?Roots.QuadraticInverse()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.QuadraticInverse()` method.


"""
quadratic_inverse(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, QuadraticInverse(); kwargs...)

superhalley(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, SuperHalley(); kwargs...)

chebyshev_like(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, ChebyshevLike(); kwargs...)

## --------------------------------------------------

## MATLAB interface to find_zero
## Main functions are
## * fzero(f, ...) to find _a_ zero of f, a univariate function
## * fzeros(f, ...) to attempt to find all zeros of f, a univariate function
## unlike `find_zero` these do not specialize on f, so
## will be faster the first use, and slower for subsequent uses (for the same f)
struct FnWrapper
    f
end
(F::FnWrapper)(x::Number) = first(F.f(x))

"""
    fzero(f, x0; order=0; kwargs...)
    fzero(f, x0, M; kwargs...)
    fzero(f, x0, M, N; kwargs...)
    fzero(f, x0; kwargs...)
    fzero(f, a::Number, b::Number; kwargs...)
    fzero(f, a::Number, b::Number; order=?, kwargs...)
    fzero(f, fp, a::Number; kwargs...)

Find zero of a function using one of several iterative algorithms.

* `f`: a scalar function or callable object

* `x0`: an initial guess, a scalar value or tuple of two values

* `order`: An integer, symbol, or string indicating the algorithm to
   use for `find_zero`. The `Order0` default may be specified directly
   by `order=0`, `order=:0`, or `order="0"`; `Order1()` by `order=1`,
   `order=:1`, `order="1"`, or `order=:secant`; `Order1B()` by
   `order="1B"`, etc.

* `M`: a specific method, as would be passed to `find_zero`, bypassing
  the use of the `order` keyword

* `N`: a specific bracketing method. When given, if a bracket is
  identified, method `N` will be used to finish instead of method `M`.

* `a`, `b`: When two values are passed along, if no `order` value is
  specified, `Bisection` will be used over the bracketing interval
  `(a,b)`. If an `order` value is specified, the value of `x0` will be set to
  `(a,b)` and the specified method will be used.

* `fp`: when `fp` is specified (assumed to compute the derivative of `f`),
  Newton's method will be used

* `kwargs...`: See `find_zero` for the specification of tolerances and other keyword arguments

Examples:
```
fzero(sin, 3)                  # use Order0() method, the default
fzero(sin, 3, order=:secant)   # use secant method (also just `order=1`)
fzero(sin, 3, Roots.Order1B()) # use secant method variant for multiple roots.
fzero(sin, 3, 4)               # use bisection method over (3,4)
fzero(sin, 3, 4, xatol=1e-6)   # use bisection method until |x_n - x_{n-1}| <= 1e-6
fzero(sin, 3, 3.1, order=1)    # use secant method with x_0=3.0, x_1 = 3.1
fzero(sin, (3, 3.1), order=2)  # use Steffensen's method with x_0=3.0, x_1 = 3.1
fzero(sin, cos, 3)             # use Newton's method
```

!!! note
    Unlike `find_zero`, `fzero` does not specialize on the type of the function argument.
    This has the advantage of making the first use of the function `f` faster, but subsequent uses slower.

"""
function fzero(f, x0::Number; kwargs...)
    x = float(x0)
    isinf(x) && throw(ConvergenceFailed("An initial value must be finite"))
    derivative_free(f, x; kwargs...)
end

function fzero(f, x0, M::AbstractUnivariateZeroMethod; kwargs...)
    find_zero(FnWrapper(f), x0, M; kwargs...)
end

function fzero(
    f,
    x0,
    M::AbstractUnivariateZeroMethod,
    N::AbstractBracketingMethod;
    kwargs...,
)
    find_zero(FnWrapper(f), x0, M, N; kwargs...)
end

function fzero(f, bracket::Tuple{T,S}; kwargs...) where {T<:Number,S<:Number}
    d = Dict(kwargs)
    if haskey(d, :order)
        find_zero(FnWrapper(f), bracket, _method_lookup[d[:order]]; kwargs...)
    else
        find_zero(FnWrapper(f), bracket, Bisection(); kwargs...)
    end
end

fzero(f, a::Number, b::Number, args...; kwargs...) = fzero(f, (a, b), args...; kwargs...)

fzero(f, x; kwargs...) = find_zero(FnWrapper(f), x; kwargs...)

fzero(f::Function, fp::Function, x0::Real; kwargs...) =
    find_zero((f, fp), x0, Newton(); kwargs...)

# match fzero up with find_zero
_method_lookup = Dict(
    0           => Order0(),
    :0          => Order0(),
    "0"         => Order0(),
    1           => Order1(),
    :1          => Order1(),
    "1"         => Order1(),
    :secant     => Order1(),
    :Secant     => Order1(),
    "1B"        => Order1B(),
    :king       => Order1B(),
    :King       => Order1B(),
    2           => Order2(),
    :2          => Order2(),
    :steffensen => Order2(),
    :Steffensen => Order2(),
    "2"         => Order2(),
    "2B"        => Order2B(),
    :esser      => Order2B(),
    :Esser      => Order2B(),
    5           => Order5(),
    :5          => Order5(),
    "5"         => Order5(),
    8           => Order8(),
    :8          => Order8(),
    "8"         => Order8(),
    16          => Order16(),
    :16         => Order16(),
    "16"        => Order16(),
)

@noinline function derivative_free(f, x0; order=0, kwargs...)
    if haskey(_method_lookup, order)
        M = _method_lookup[order]
    else
        warn("Invalid order specified. See ?fzero.")
        throw(ArgumentError())
    end

    # d = (kv[1] == :ftol ? :atol=>kv[2] :
    #      kv[1] == :ftolrel ? :rtol=>kv[2] :
    #      kv[1] == :xtol ? :xatol=>kv[2] :
    #      kv[1] == :xtolrel ? xrtol=>kv[2] :
    #      kv[1] => kv[1] for kv in kwargs)

    d = Dict(kwargs)
    for (o, n) in ((:ftol, :atol), (:ftolrel, :rtol), (:xtol, :xatol), (:xtolrel, :xrtol))
        if haskey(d, o)
            d[n] = d[o]
        end
    end

    find_zero(FnWrapper(f), x0, M; d...)
end

## fzeros
"""

`fzeros(f, a, b; kwargs...)`

Searches for all zeros of `f` within an interval `(a,b)`. Assume neither `a` or `b` is a zero.

Dispatches to `find_zeros(f, a, b; kwargs...)`.
"""
function fzeros(f, a::Number, b::Number; kwargs...)
    find_zeros(FnWrapper(f), float(a), float(b); kwargs...)
end
fzeros(f, bracket::Vector{T}; kwargs...) where {T<:Number} =
    fzeros(f, bracket[1], bracket[2]; kwargs...)
fzeros(f, bracket::Tuple{T,S}; kwargs...) where {T<:Number,S<:Number} =
    fzeros(f, bracket[1], bracket[2]; kwargs...)
