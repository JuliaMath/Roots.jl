##################################################
## Classical derivative-based, iterative, root-finding algorithms
##
## If ri = f^(i-1)/f^(i), then these have an update step `xn - delta` where:
##
## * Newton: delta = r1  # order 2 if simple root (multiplicity 1)
## * Halley: delta = 2*r2/(2r2 - r1) * r1 # order 3 if simple root
## * Schroder: delta = r2  / (r2 - r1) * r1  # order 2
## * Thukral(3): delta =  (-2*r3)*(r2 - r1)/(r1^2 - 3*r1*r3 + 2*r2*r3) * r1 # order 3
## * Thukral(4): delta =  3*r1*r2*r4*(r1^2 - 3*r1*r3 + 2*r2*r3)/(-r1^3*r2 + 4*r1^2*r2*r4 + 3*r1^2*r3*r4 - 12*r1*r2*r3*r4 + 6*r2^2*r3*r4) # order 4
##
## The latter two come from
## [Thukral](http://article.sapub.org/10.5923.j.ajcam.20170702.05.html). They are not implemented.

## Newton
abstract type AbstractNewtonLikeMethod <: AbstractUnivariateZeroMethod end
fn_argout(::AbstractNewtonLikeMethod) = 2
struct Newton <: AbstractNewtonLikeMethod end
"""

    Roots.Newton()

Implements Newton's [method](http://tinyurl.com/b4d7vls):
`xᵢ₊₁ =  xᵢ - f(xᵢ)/f'(xᵢ)`.  This is a quadratically convergent method requiring
one derivative. Two function calls per step.

Example

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin,cos), 3.0, Roots.Newton()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which returns (f, f/f') as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x)), 3.0, Roots.Newton()) ≈ π
true
```

This can be advantageous if the derivative is easily computed from the
value of f, but otherwise would be expensive to compute.

The error, `eᵢ = xᵢ - α`, can be expressed as `eᵢ₊₁ = f[xᵢ,xᵢ,α]/(2f[xᵢ,xᵢ])eᵢ²` (Sidi, Unified treatment of regula falsi, Newton-Raphson, secant, and Steffensen methods for nonlinear equations).

"""
Newton

# we store x0,x1,fx0,fx1 **and** Δ = fx1/f'(x1)
struct NewtonState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δ::T
    fxn1::S
    fxn0::S
end

function init_state(M::Newton, F::Callable_Function, x)
    x₀ = float(first(x))
    T = eltype(x₀)
    fx₀, Δ::T = F(x₀)
    x₁::T = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

# compute fx₁, Δ
function init_state(::Newton, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    fx₁, Δ::T = F(x₁)
    NewtonState(x₁, x₀, Δ, fx₁, fx₀)
end

initial_fncalls(M::Newton) = 2

function update_state(
    method::Newton,
    F,
    o::NewtonState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ::T = o.Δ

    if isissue(Δ)
        log_message(l, "Issue with `f/f′'")
        return o, true
    end

    xn0, xn1 = xn1, xn1 - Δ
    fxn0 = fxn1
    fxn1, Δ = F(xn1)
    incfn(l, 2)

    @set! o.xn0 = xn0
    @set! o.xn1 = xn1
    @set! o.Δ = Δ
    @set! o.fxn0 = fxn0
    @set! o.fxn1 = fxn1

    return o, false
end

"""
    Roots.newton(f, fp, x0; kwargs...)

Implementation of Newton's method: `xᵢ₊₁ =  xᵢ - f(xᵢ)/f'(xᵢ)`.

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`.

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero` using the `Roots.Newton()` method.

See also `Roots.newton((f,fp), x0)` and `Roots.newton(fΔf, x0)` for simpler implementations.

"""
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)

## Halley
abstract type AbstractHalleyLikeMethod <: AbstractUnivariateZeroMethod end
abstract type AbstractΔMethod <: AbstractHalleyLikeMethod end

fn_argout(::AbstractHalleyLikeMethod) = 3

"""
    Roots.Halley()

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 - (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)^(-1)`
This method is cubically converging, but requires more function calls per step (3) than
other methods.

Example

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Halley()) ≈ π
true
```

This can be advantageous if the derivatives are easily computed from
the computation for f, but otherwise would be expensive to compute separately.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ ≈ -(2f'⋅f''' -3⋅(f'')²)/(12⋅(f'')²) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct Halley <: AbstractΔMethod end

struct HalleyState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δ::T
    ΔΔ::T
    fxn1::S
    fxn0::S
end

# we compute one step here to get x₁
function init_state(M::AbstractHalleyLikeMethod, F::Callable_Function, x)
    x₀ = float(first(x))
    T = eltype(x₀)
    fx₀, (Δ, ΔΔ) = F(x₀)
    Δx = calculateΔ(M, Δ, ΔΔ)
    x₁::T = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

function init_state(::AbstractHalleyLikeMethod, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    fx₁, (Δ::T, ΔΔ::T) = F(x₁)
    HalleyState(x₁, x₀, Δ, ΔΔ, fx₁, fx₀)
end

initial_fncalls(M::AbstractHalleyLikeMethod) = 2 * 3

calculateΔ(method::Halley, r1, r2) = 2r2 / (2r2 - r1) * r1

function update_state(
    method::AbstractΔMethod,
    F,
    o::HalleyState{T,S},
    options::UnivariateZeroOptions,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ::T = calculateΔ(method, r1, r2)
    if isissue(Δ)
        log_message(l, "Issue with computing `Δ`")
        return (o, true)
    end

    xn1::T = xn - Δ
    fxn1::S, (r1::T, r2::T) = F(xn1)
    incfn(l, 3)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1
    @set! o.Δ = r1
    @set! o.ΔΔ = r2

    return o, false
end

"""
    Roots.halley(f, fp, fpp, x0; kwargs...)

Implementation of Halley's method (cf `?Roots.Halley()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.Halley()` method.


"""
halley(f, fp, fpp, x0; kwargs...) = find_zero((f, fp, fpp), x0, Halley(); kwargs...)

"""
    Roots.QuadraticInverse()

Implements the [quadratic inverse method](https://doi.org/10.2307/2322644) also known as [Chebyshev's method]((https://dl.acm.org/doi/10.1080/00207160802208358)),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 + (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)`.
This method is cubically converging, but requires more function calls per step (3) than
other methods.

Example

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.QuadraticInverse()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.QuadraticInverse()) ≈ π
true
```

This can be advantageous if the derivatives are easily computed from
the computation for f, but otherwise would be expensive to compute separately.

The error, `eᵢ = xᵢ - α`, [satisfies](https://dl.acm.org/doi/10.1080/00207160802208358)
`eᵢ₊₁ ≈ (1/2⋅(f''/f')² - 1/6⋅f'''/f')) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct QuadraticInverse <: AbstractΔMethod end

calculateΔ(method::QuadraticInverse, r1, r2) = (1 + r1 / (2r2)) * r1

"""
    Roots.quadratic_inverse(f, fp, fpp, x0; kwargs...)

Implementation of the quadratic inverse method (cf `?Roots.QuadraticInverse()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.QuadraticInverse()` method.


"""
quadratic_inverse(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, QuadraticInverse(); kwargs...)

"""
CHEBYSHEV-LIKE METHODS AND QUADRATIC EQUATIONS (J. A. EZQUERRO, J. M. GUTIÉRREZ, M. A. HERNÁNDEZ and M. A. SALANOVA)
"""
struct ChebyshevLike <: AbstractΔMethod end

calculateΔ(method::ChebyshevLike, r1, r2) = (1 + r1 / (2r2) * (1 + r1 / r2)) * r1

chebyshev_like(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, ChebyshevLike(); kwargs...)

"""
An acceleration of Newton's method: Super-Halley method (J.M. Gutierrez, M.A. Hernandez
"""
struct SuperHalley <: AbstractΔMethod end

calculateΔ(method::SuperHalley, r1, r2) = (1 + r1 / (2r2 - 2r1)) * r1

superhalley(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, SuperHalley(); kwargs...)
