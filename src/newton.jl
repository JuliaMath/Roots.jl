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
"""

    Roots.Newton()

Implements Newton's [method](http://tinyurl.com/b4d7vls): `x_n1 = xn -
f(xn)/f'(xn)`.  This is a quadratically converging method requiring
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

The error, `en = xn - alpha`, can be expressed as `e1 = f[x0,x0,alpha]/(2f[x0,x0])e0^2` (Sidi).

"""
struct Newton <: AbstractNewtonLikeMethod end

function init_state(method::AbstractNewtonLikeMethod, fs, x)

    x1 = float(x)
    T = eltype(x1)
    tmp = fΔx(fs, x1)
    fx1, Δ::T = tmp[1], tmp[2] # faster to pass in fx, fx/f'(x) than (fx, f'(x)) and compute

    fnevals = 1
    S = eltype(fx1)

    state = UnivariateZeroState(x1, oneunit(x1) * (0*x1)/(0*x1), [Δ],
                                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[],
                                0, fnevals,
                                false, false, false, false,
                                "")
    state
end

function init_state!(state::UnivariateZeroState{T,S}, M::Newton, fs, x) where {T,S}
    x1::T = float(x)
    tmp = fΔx(fs, x)
    fx1::S, Δ::T = tmp[1], tmp[2]

     init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), [Δ],
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end


function update_state(method::Newton, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    xn = o.xn1
    fxn = o.fxn1
    r1 = o.m[1]

    if isissue(r1)
        o.stopped=true
        return
    end

    xn1 = xn - r1

    tmp = fΔx(fs, xn1)
    fxn1::S, r1::T = tmp[1], tmp[2]
    incfn(o,2)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    empty!(o.m); push!(o.m, r1)

    nothing



end

"""
    Roots.newton(f, fp, x0; kwargs...)

Implementation of Newton's method: `x_n1 = x_n - f(x_n)/ f'(x_n)`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`.

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero` using the `Roots.Newton()` method.

See also `Roots.newton((f,fp), x0) and `Roots.newton(fΔf, x0)` for simpler implementations.

"""
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)


## Halley
abstract type AbstractHalleyLikeMethod <: AbstractUnivariateZeroMethod end

"""
    Roots.Halley()

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`x_n1 = xn - f/f' * (1 - f/f' * f''/f' * 1/2)^(-1)
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

The error, `e_n = x_n - alpha`, satisfies
`e1 ≈ -(2f'⋅f''' -3⋅(f'')^2)/(12⋅(f'')^2) ⋅ e0^3` (all evaluated at `alpha`).

"""
struct Halley <: AbstractHalleyLikeMethod
end


function init_state(method::AbstractHalleyLikeMethod, fs, x)

    x1 = float(x)
    T = eltype(x1)
    tmp = fΔxΔΔx(fs, x1)
    fx1, Δ::T, ΔΔ::T = tmp[1], tmp[2], tmp[3]
    S = eltype(fx1)
    fnevals = 3

    state = UnivariateZeroState(x1, oneunit(x1) * (0*x1)/(0*x1), [Δ,ΔΔ],
                                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[],
                                0, fnevals,
                                false, false, false, false,
                                "")
    state
end

function init_state!(state::UnivariateZeroState{T,S}, M::AbstractHalleyLikeMethod, fs, x) where {T,S}
    x1::T = float(x)
    tmp = fΔxΔΔx(fs, x)
    fx1::S, Δ::T, ΔΔ::T = tmp[1], tmp[2], tmp[3]

     init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), [Δ, ΔΔ],
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end

function update_state(method::Halley, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.m

    xn1::T = xn - 2*r2/(2r2 - r1) * r1

    tmp = fΔxΔΔx(fs, xn1)
    fxn1::S, r1::T, r2::T = tmp[1], tmp[2], tmp[3]
    incfn(o,3)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    empty!(o.m); append!(o.m, (r1, r2))
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
#halley(f,  x0; kwargs...) = find_zero(f, x0, Halley(); kwargs...) # deprecated
#halley(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Halley(); kwargs...) # deprecated
@deprecate halley(f,  x0; kwargs...)    halley(f, fp, fpp, x0; kwargs...)
@deprecate halley(f, fp, x0; kwargs...) halley(f, fp, fpp, x0; kwargs...)


## Shroder-like methods

"""
    Roots.Schroder()


Schröder's method, like Halley's method, utilizes f, f', and
f''. Unlike Halley it is quadratically converging, but this is
independent of the multiplicity of the zero (cf. Schröder, E. "Über
unendlich viele Algorithmen zur Auflösung der Gleichungen."
Math. Ann. 2, 317-365, 1870;
http://mathworld.wolfram.com/SchroedersMethod.html). (Schröder's
method applies Newton's method to `f/f'`, a function with all
simple zeros.)


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
returns (f, f/f',f'/f'') as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Schroder())
```

This can be advantageous if the derivatives are easily computed from
the value of f, but otherwise would be expensive to compute.

The error, `e_n = x_n - alpha`, is the same as `Newton` with `f` replaced by `f/f'`.

"""
struct Schroder <: AbstractHalleyLikeMethod
end
const Schroeder = Schroder # either spelling
const Schröder = Schroder

function update_state(method::Schroder, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.m[1], o.m[2]

    delta =  r2 / (r2 - r1) * r1  # m * r1
    xn1::T = xn - delta

    tmp = fΔxΔΔx(fs, xn1)
    fxn1::S, r1::T, r2::T = tmp[1], tmp[2], tmp[3]
    incfn(o,3)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    o.m[1], o.m[2] = r1, r2
end
