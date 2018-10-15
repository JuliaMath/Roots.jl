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
    #fx1, Δ::T = fΔx(fs, x1)
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
struct Halley <: AbstractUnivariateZeroMethod
end


function init_state(method::Halley, fs, x)

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

function init_state!(state::UnivariateZeroState{T,S}, M::Halley, fs, x) where {T,S}
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
