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
```
find_zero((sin,cos), 3.0, Roots.Newton())
```

If function evaluations are expensive one can pass in a function which returns (f, f/f') as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x)), 3.0, Roots.Newton())
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
    fx₀, Δ = F(x₀)
    x₁ = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

# compute fx₁, Δ
function init_state(::Newton, F, x₀, x₁, fx₀, fx₁)
    fx₁, Δ = F(x₁)
    NewtonState(x₁, x₀, Δ, fx₁, fx₀)
end

initial_fncalls(M::Newton) = 2


function update_state(method::Newton, F, o::NewtonState{T,S}, options, l=NullTracks()) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ = o.Δ

    if isissue(Δ)
        log_message(l, "Issue with `f/f′'")
        return o, true
    end

    xn0, xn1 = xn1, xn1-Δ
    fxn0 = fxn1
    fxn1, Δ = F(xn1)
    incfn(l,2)

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
fn_argout(::AbstractHalleyLikeMethod) = 3

"""
    Roots.Halley()

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 - (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)^(-1)`
This method is cubically converging, but requires more function calls per step (3) than
other methods.

Example
```
find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley())
```

If function evaluations are expensive one can pass in a function which
returns (f, f/f',f'/f'') as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Halley())
```

This can be advantageous if the derivatives are easily computed from
the value of f, but otherwise would be expensive to compute.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ ≈ -(2f'⋅f''' -3⋅(f'')²)/(12⋅(f'')²) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct Halley <: AbstractHalleyLikeMethod
end

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
    fx₀, Δ, ΔΔ = F(x₀)
    x₁ = x₀ - 2ΔΔ/(2ΔΔ - Δ) * Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

function init_state(::AbstractHalleyLikeMethod, F, x₀, x₁, fx₀, fx₁)
    fx₁, Δ, ΔΔ = F(x₁)
    HalleyState(x₁, x₀, Δ, ΔΔ, fx₁, fx₀)
end

initial_fncalls(M::AbstractHalleyLikeMethod) = 2*3

function update_state(method::Halley, F, o::HalleyState{T,S}, options::UnivariateZeroOptions, l=NullTracks()) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ =  2*r2/(2r2 - r1) * r1
    if isissue(Δ)
        log_message(l, "Issue with computing `Δ`")
        return (o, true)
    end

    xn1::T = xn - Δ
    fxn1::S, r1::T, r2::T = F(xn1)
    incfn(l,3)


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
    Roots.Schroder()


Schröder's method, like Halley's method, utilizes `f`, `f'`, and
`f''`. Unlike Halley it is quadratically converging, but this is
independent of the multiplicity of the zero (cf. Schröder, E. "Über
unendlich viele Algorithmen zur Auflösung der Gleichungen."
Math. Ann. 2, 317-365, 1870;
[mathworld](http://mathworld.wolfram.com/SchroedersMethod.html)). Schröder's
method applies Newton's method to `f/f'`, a function with all
simple zeros.


Example
```
m = 2
f(x) = (cos(x)-x)^m
fp(x) = (-x + cos(x))*(-2*sin(x) - 2)
fpp(x) = 2*((x - cos(x))*cos(x) + (sin(x) + 1)^2)
find_zero((f, fp, fpp), pi/4, Roots.Halley())    # 14 steps
find_zero((f, fp, fpp), 1.0, Roots.Schroder())    # 3 steps
```

(Whereas, when `m=1`, Halley is 2 steps to Schröder's 3.)

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Schroder())
```

This can be advantageous if the derivatives are easily computed from
the value of `f`, but otherwise would be expensive to compute.

The error, `eᵢ = xᵢ - α`, is the same as `Newton` with `f` replaced by `f/f'`.

"""
struct Schroder <: AbstractHalleyLikeMethod
end
const Schroeder = Schroder # either spelling
const Schröder = Schroder

## Shroder-like methods
function init_state(M::Schroder, F::Callable_Function, x)
    x₀ = float(first(x))
    fx₀, Δ, ΔΔ = F(x₀)
    x₁ = x₀ - ΔΔ/(ΔΔ - Δ) * Δ # m*r1
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end



function update_state(method::Schroder, F, o::AbstractUnivariateZeroState{T,S},
                      options::UnivariateZeroOptions, l=NullTracks()) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ =  r2 / (r2 - r1) * r1  # m * r1

    if isissue(Δ)
        log_message(l, "Issue with increment")
        return o, true
    end

    xn1::T = xn - Δ

    fxn1::S, r1::T, r2::T = F(xn1)
    incfn(l,3)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1
    @set! o.Δ = r1
    @set! o.ΔΔ = r2

    return o, false
end
