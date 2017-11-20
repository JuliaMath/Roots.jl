##################################################
## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley
## Historic, we have derivative free versions of similar order


## Newton
"""

    Roots.Newton()

Implements Newton's [method](http://tinyurl.com/b4d7vls): `x_n1 = xn -
f(xn)/f'(xn)`.  This is a quadratically converging method requiring
one derivative. If a derivative is not specified, the `ForwardDiff` package
will be used, as applicable.

Unlike other methods, this method accepts complex inputs.
"""
struct Newton <: UnivariateZeroMethod
end

function callable_function(method::Newton, f::Tuple, x0)
    length(f) == 1 && return FirstDerivative(f[1], D(f[1]))
    FirstDerivative(f[1], f[2], f[1](x0))
end
callable_function(method::Newton, f::Any, x0) = FirstDerivative(f, D(f), x0)

function update_state(method::Newton, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}
    xn = o.xn1
    fxn = o.fxn1
    fpxn = fs.fp(xn)

    if isissue(fpxn)
        o.stopped=true
        return
    end
    
    xn1 = xn - fxn / fpxn
    fxn1 = fs.f(xn1)
    incfn(o)
    
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1


    incsteps(o)
    
end

## extra work to allow for complex values

## extra work to allow for complex values
function derivative_free_setup(method::Newton, fs::CallableFunction, x0::T;
                                  bracket=missing,
                                  xabstol=zero(T), xreltol=zero(T),
                                  abstol=4*eps(T), reltol=4*eps(T),
                                  maxevals=40, maxfnevals=typemax(Int),
                                  verbose::Bool=false) where {T<:AbstractFloat}
    x = float(x0)


    prob = UnivariateZeroProblem(fs, x, bracket)
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol,  maxevals, maxfnevals, verbose)
    prob, options
end

function derivative_free_setup(method::Newton, fs::CallableFunction, x0::Complex{T};
                                  bracket=missing,
                                  xabstol=zero(T), xreltol=zero(T),
                                  abstol=4*eps(T), reltol=4*eps(T),
                                  maxevals=40, maxfnevals=typemax(Int),
                                  verbose::Bool=false) where {T<:AbstractFloat}
    x = float(x0)
    bracket = missing     # bracket makes no sense for complex input, but one is expected

    prob = UnivariateZeroProblem(fs, x, bracket)
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol,  maxevals, maxfnevals, verbose)
    prob, options
end

"""

Implementation of Newton's method: `x_n1 = x_n - f(x_n)/ f'(x_n)`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function=D(f)` -- derivative of `f`. Defaults to automatic derivative

* `x0::Number` -- initial guess. For Newton's method this may be complex.

Keyword arguments:

* `ftol`. Stop iterating when |f(xn)| <= max(1, |xn|) * ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `maxfneval`. Stop iterating if more than this many function calls, throw error.

* `verbose::Bool=false` Set to `true` to see trace.

"""
newton(f, x0; kwargs...) = find_zero(f, x0, Newton(); kwargs...)
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)


## Halley


"""
    
    Roots.Halley()

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`x_n1 = xn - (2 f(xn)*f'(xn)) / (2 f'(xn)^2 - f(xn) * f''(xn))`.
This method is cubically converging, but requires more function calls per step than
other methods.
"""    
struct Halley <: UnivariateZeroMethod
end

function callable_function(method::Halley, f::Tuple, x0)
    length(f) == 1 && return SecondDerivative(f[1], D(f[1]), D(f[1],2))
    length(f) == 2 && return SecondDerivative(f[1], f[2], D(f[2],1))
    SecondDerivative(f[1], f[2], f[3], f[1](x0))
end
callable_function(method::Halley, f, x0) = SecondDerivative(f, D(f), D(f, 2))


function update_state(method::Halley, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}
    xn = o.xn1
    fxn = o.fxn1
    fpxn = fs.fp(xn); incfn(o)
    fppxn = fs.fpp(xn); incfn(o)
    
    xn1 = xn - 2fxn*fpxn / (2*fpxn*fpxn - fxn * fppxn)
    fxn1 = fs.f(xn1); incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    incsteps(o)
end

"""

Implementation of Halley's method. `xn1 = xn - 2f(xn)*f'(xn) / (2*f'(xn)^2 - f(xn) * f''(xn))`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function=D(f)` -- derivative of `f`. Defaults to automatic derivative

* `fpp:Function=D(f,2)` -- second derivative of `f`.

* `x0::Real` -- initial guess

Keyword arguments:

* `ftol`. Stop iterating when |f(xn)| <= max(1, |xn|) * ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `verbose::Bool=false` Set to `true` to see trace.

"""
halley(f,  x0; kwargs...) = find_zero(f, x0, Halley(); kwargs...)
halley(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Halley(); kwargs...)
halley(f, fp, fpp, x0; kwargs...) = find_zero((f, fp, fpp), x0, Halley(); kwargs...)
