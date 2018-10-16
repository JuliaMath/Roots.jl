##################################################
## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley
## Historic, we have derivative free versions of similar order


## Newton
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

"""
struct Newton <: AbstractUnivariateZeroMethod
end

function init_state(method::Newton, fs, x)

    x1 = float(x)
    T = eltype(x1)
    tmp = fΔx(fs, x1)
    fx1, Δ::T = tmp[1], tmp[2]

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
    fx1::S, Δ::T = tmp[1],tmp[2]

     init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), [Δ],
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end


function update_state(method::Newton, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    xn = o.xn1
    fxn = o.fxn1
    Δxn::T = o.m[1]

    if isissue(Δxn)
        o.stopped=true
        return
    end

    xn1 = xn - Δxn

    fxn1::S, Δxn1::T = fΔx(fs, xn1)
    incfn(o,2)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    empty!(o.m); push!(o.m, Δxn1)

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
    fx1::S, Δ::T, ΔΔ::T = tmp[1],tmp[2], tmp[3]

     init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), [Δ],
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end

function update_state(method::Halley, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    Δxn, ΔΔxn = o.m


    xn1::T = xn - Δxn * inv(1  - Δxn / ΔΔxn / 2)

    fxn1::S, Δxn1::T, ΔΔxn1::T = fΔxΔΔx(fs, xn1)
    incfn(o,3)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    empty!(o.m); append!(o.m, (Δxn1, ΔΔxn1))
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

## a method needing derivatives
abstract type AbstractSchroderLikeMethod <: AbstractUnivariateZeroMethod end

#
fΔxs(F::DerivativeFree, x) = F.f(x)
function fΔxs(F::FirstDerivative, x)
    u, v = F.f(x), F.fp(x)
    (u, u/v)
end
function fΔxs(F::SecondDerivative, x)
    u, v, w = F.f(x), F.fp(x), F.fpp(x)
    (u, u/v, v/w)
end
fΔxs(F, x) = F(x)



peel(xs, T) = (xs[1], collect(T, xs[2:end]))
peel(xs) = (xs[1], collect(xs[2:end]))

# can speed these up if they are specific like above cases
function init_state(method::AbstractSchroderLikeMethod, fs, x)

    x1 = float(x)
    T = eltype(x1)
    fx1, Δs = peel(fΔxs(fs, x1), T)
    S = eltype(fx1)
    fnevals = 3


    state = UnivariateZeroState(x1, oneunit(x1) * (0*x1)/(0*x1), Δs,
                                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[],
                                0, fnevals,
                                false, false, false, false,
                                "")
    state
end

function init_state!(state::UnivariateZeroState{T,S}, M::AbstractSchroderLikeMethod, fs, x) where {T,S}
    x1::T = float(x)
    fx1::S, Δs = peel(fΔxs(fs, x1), T)


    init_state!(state, x1, oneunit(x1) * (0*x1)/(0*x1), Δs,
                fx1, oneunit(fx1) * (0*fx1)/(0*fx1), S[])
end


"""
    Roots.Schroder()


Schroder's method, like Halley's method, utilizes f, f', and
f''. Unlike Halley it is quadratically converging, but this is independent of the
multiplicity of the zero.

If ri = f^(i-1)/f^(i), then we have an update step `xn - delta` where:

* Newton: delta = r1
* Halley: delta = 2*r2*r2/(2r2 - r1)
* Schroder: delta = r2 * r1 / (r2 - r1)
* Thukral3: delta =  -2*r1*r3*(r2 - r1)/(r1^2 - 3*r1*r3 + 2*r2*r3)
* Thukral4: delta =  3*r1*r2*r4*(r1^2 - 3*r1*r3 + 2*r2*r3)/(-r1^3*r2 + 4*r1^2*r2*r4 + 3*r1^2*r3*r4 - 12*r1*r2*r3*r4 + 6*r2^2*r3*r4)

The latter two come from
[Thukral](http://article.sapub.org/10.5923.j.ajcam.20170702.05.html)


Example
```
m = 2
f(x) = (cos(x)-x)^m
fp(x) = (-x + cos(x))*(-2*sin(x) - 2)
fpp(x) = 2*((x - cos(x))*cos(x) + (sin(x) + 1)^2)
find_zero((f, fp, fpp), pi/4, Roots.Halley())    # 14 steps
find_zero((f, fp, fpp), 1.0, Roots.Schroder())    # 3 steps
```

(Whereas, when `m=1`, Halley is 2 steps to Schroder's 3.)

If function evaluations are expensive one can pass in a function which
returns (f, f/f',f'/f'') as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Schroder())
```

This can be advantageous if the derivatives are easily computed from
the value of f, but otherwise would be expensive to compute.

"""
struct Schroder <: AbstractHalleyLikeMethod
end


function update_state(method::Schroder, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.m[1], o.m[2]

    delta =  r2 / (r2 - r1) * r1  # m * r1
    xn1::T = xn - delta

    fxn1::S, r1::T, r2::T = fΔxΔΔx(fs, xn1)
    incfn(o,3)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    empty!(o.m); append!(o.m, (r1, r2))
end


"""
    Roots.Thukral3()

Higher order Schroder-like algorithm. Requires 3 derivatives. Trades
off one order of convergence for an algorithm that does not depend on
the multiplicity of the zero.

From [Thukral](http://article.sapub.org/10.5923.j.ajcam.20170702.05.html).


Example
```
using ForwardDiff
D(f, n=1) = n > 1 ? D(D(f),n-1) : x -> ForwardDiff.derivative(f, float(x))

m = 5
f(x) = (exp(-x^2) - exp(x^2) - x^8 + 10)^m
function Fs(f, x)
    fx = f(x)
    fp, fpp, fppp, fpppp = [D(f,i)(x) for i in 1:4]
    r1, r2, r3, r4 = fx/fp, fp/fpp, fpp/fppp, fppp/fpppp
    fx, r1, r2, r3, r4
end

x0 = 1.4

find_zero(x->Fs(f, x), x0, Roots.Thukral4()) # 2 iterations
find_zero(x->Fs(f, x), x0, Roots.Thukral3()) # 2 iterations
find_zero(x->Fs(f, x), x0, Roots.Halley())   # 24 iterations
```
"""
struct Thukral3 <: AbstractSchroderLikeMethod
end


function update_state(method::Thukral3, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1::T, r2::T, r3::T = o.m[1], o.m[2], o.m[3]

    delta::T = -2*r1*r3*(r1 - r2)/(r1^2 - 3*r1*r3 + 2*r2*r3)

    xn1::T = xn - delta


    fxn1::S, Δs::Vector{T} = peel(fΔxs(fs, xn1), T)
    r1,r2,r3 = Δs[1], Δs[2], Δs[3]
    incfn(o,4)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    o.m[1], o.m[2], o.m[3] = r1, r2, r3

end

"""
    Roots.Thukral4()

Higher order Schroder-like algorithm. Requires 4 derivatives. Trades
off one order of convergence for an algorithm that does not depend on
the multiplicity of the zero.

From [Thukral](http://article.sapub.org/10.5923.j.ajcam.20170702.05.html).


Example
```
using ForwardDiff
D(f, n=1) = n > 1 ? D(D(f),n-1) : x -> ForwardDiff.derivative(f, float(x))

m = 5
f(x) = (exp(-x^2) - exp(x^2) - x^8 + 10)^m
function Fs(f, x)
    fx = f(x)
    fp, fpp, fppp, fpppp = [D(f,i)(x) for i in 1:4]
    r1, r2, r3, r4 = fx/fp, fp/fpp, fpp/fppp, fppp/fpppp
    fx, r1, r2, r3, r4
end

x0 = 1.4

find_zero(x->Fs(f, x), x0, Roots.Thukral4()) # 2 iterations
find_zero(x->Fs(f, x), x0, Roots.Thukral3()) # 2 iterations
find_zero(x->Fs(f, x), x0, Roots.Halley())   # 24 iterations
```
"""
struct Thukral4 <: AbstractSchroderLikeMethod
end


function update_state(method::Thukral4, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2, r3, r4 = o.m[1], o.m[2], o.m[3], o.m[4]

    delta::T = 3*r1*r2*r4*(r1^2 - 3*r1*r3 + 2*r2*r3)/(-r1^3*r2 + 4*r1^2*r2*r4 + 3*r1^2*r3*r4 - 12*r1*r2*r3*r4 + 6*r2^2*r3*r4)
    xn1::T = xn - delta

    fxn1::S, Δs::Vector{T} = peel(fΔxs(fs, xn1), T)
    r1,r2,r3, r4 = Δs[1], Δs[2], Δs[3], Δs[4]
    incfn(o,5)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    o.m[1], o.m[2], o.m[3], o.m[4] = r1, r2, r3, r4
end
