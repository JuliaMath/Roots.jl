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
find_zero(Roots.fg(x -> (sin(x), sin(x)/cos(x))), 3.0, Roots.Newton())
```

This can be advantageous if the derivative is easily computed from the
value of f, but otherwise would be expensive to compute.

"""
struct Newton <: AbstractUnivariateZeroMethod
end

function update_state(method::Newton, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    xn = o.xn1
    #fxn = o.fxn1
    #fpxn = fs(xn,1)
    fxn::T, Δxn::S = fΔf(fs, xn)

    if isissue(Δxn)
        o.stopped=true
        return
    end

    xn1 = xn - Δxn
    fxn1 = fxn#fs(xn1)
    incfn(o)
    incfn(o)
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1




end

"""
    newton(f, fp, x0; kwargs...)

Implementation of Newton's method: `x_n1 = x_n - f(x_n)/ f'(x_n)`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`.

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero` using the `Roots.Newton()` method.

See also `newton((f,fp), x0) and `newton(fΔf, x0)` for simpler implementations.

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

If function evaluations are expensive one can pass in a function which returns (f, f/f',f''/f) as follows

```
find_zero(Roots.fg(x -> (sin(x), sin(x)/cos(x), -sin(x)/cos(x))), 3.0, Roots.Halley())
```

This can be advantageous if the derivatives are easily computed from
the value of f, but otherwise would be expensive to compute.

"""
struct Halley <: AbstractUnivariateZeroMethod
end


function update_state(method::Halley, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
#    fxn = o.fxn1
#    fpxn = fs(xn,1); incfn(o)
    #    fppxn = fs(xn,2); incfn(o)
    fxn::T, Δxn::S, ΔΔxni::S = fΔfΔΔf(fs, xn)  # (f, f/f', f''/f')
    incfn(o,3)

    xn1 = xn - Δxn * inv(1  - Δxn * ΔΔxni/2)

#    fxn1 = fs(xn1); incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn # lag one
end

"""
    halley(f, fp, fpp, x0; kwargs...)

Implementation of Halley's method. `xn1 = xn - 2f(xn)*f'(xn) / (2*f'(xn)^2 - f(xn) * f''(xn))`

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
