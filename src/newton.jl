##################################################
## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley
## Historic, we have derivative free versions of similar order


## Newton
"""

    Roots.Newton()

Implements Newton's [method](http://tinyurl.com/b4d7vls): `x_n1 = xn -
f(xn)/f'(xn)`.  This is a quadratically converging method requiring
one derivative. A derivative need not be specified, as the `ForwardDiff` package
may be used to compute this.
"""
struct Newton <: AbstractUnivariateZeroMethod
end

function update_state(method::Newton, fs, o, options) 
    xn = o.xn1
    fxn = o.fxn1
    fpxn = fs(xn,1)

    if isissue(fpxn)
        o.stopped=true
        return
    end
    
    xn1 = xn - fxn / fpxn
    fxn1 = fs(xn1)
    incfn(o)
    
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1


    incsteps(o)
    
end

"""

Implementation of Newton's method: `x_n1 = x_n - f(x_n)/ f'(x_n)`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`. 

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, x)` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero`.

"""
newton(f, x0; kwargs...) = find_zero(f, x0, Newton(); kwargs...) # deprecated now
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)


## Halley


"""
    
    Roots.Halley()

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`x_n1 = xn - (2 f(xn)*f'(xn)) / (2 f'(xn)^2 - f(xn) * f''(xn))`.
This method is cubically converging, but requires more function calls per step than
other methods.
"""    
struct Halley <: AbstractUnivariateZeroMethod
end


function update_state(method::Halley, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    fpxn = fs(xn,1); incfn(o)
    fppxn = fs(xn,2); incfn(o)
    
    xn1 = xn - 2fxn*fpxn / (2*fpxn*fpxn - fxn * fppxn)
    fxn1 = fs(xn1); incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1
    incsteps(o)
end

"""

Implementation of Halley's method. `xn1 = xn - 2f(xn)*f'(xn) / (2*f'(xn)^2 - f(xn) * f''(xn))`
    
Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, x)` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero`.

"""
halley(f,  x0; kwargs...) = find_zero(f, x0, Halley(); kwargs...)
halley(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Halley(); kwargs...)
halley(f, fp, fpp, x0; kwargs...) = find_zero((f, fp, fpp), x0, Halley(); kwargs...)
