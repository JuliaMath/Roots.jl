## A collection of historical methods for pedagogical purposes.
##
## secant_method, newton, and halley
##
## These have an argument `verbose` that can be specified to get a trace of the algorithm

function fn_template(meth, args...;
                     verbose::Bool=false,kwargs...)
    
    out = meth(args...; verbose=verbose, kwargs...)
    
    if verbose
        verbose_output(out)
    end
    
    if out.state == :converged
        out.x[end]
    else
        steps = out.iterations > 1 ? "steps" : "step"
        throw(ConvergenceFailed("Failed to converge in $(out.iterations) $steps."))
    end
end

## order 1 secant method

"""

Implementation of secant method: `x_n1 = x_n - f(x_n) * f(x_n)/ (f(x_n) - f(x_{n-1}))`

Arguments:

* `f::Function` -- function to find zero of

* `x0::Real` -- initial guess is [x0, x1]

* `x1::Real` -- initial guess is [x0, x1]

Keyword arguments:

* `ftol`. Stop iterating when |f(xn)| <= max(1, |xn|) * ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `maxfneval`. Stop iterating if more than this many function calls, throw error.

* `verbose::Bool=false` Set to `true` to see trace.

"""
secant_method(f::Function, x0::Real, x1::Real;
              kwargs...) = fn_template(secmeth, f, x0, x1; kwargs...)




## Add derivatives for newton, halley
type ZeroFunction3{S<:Number, T<:FloatingPoint} <: ZeroFunction
    f
    fp
    fpp
    x::Vector{S}
    fxn::S
    update::Function
    state::Symbol
    how::Symbol
    iterations::Int
    fncalls::Int
    ftol::T
    xtol::T
    xtolrel::T
end


# Newton-Raphson method (quadratic convergence)
function newton_update(F)
    x0 = F.x[end]
    x1 = x0 - F.fxn / F.fp(x0)

    F.fncalls += 2
    F.fxn = F.f(x1)
    push!(F.x, x1)
end
 

# this version uses ZeroFunctionComp for Complex numbers
function newtonmeth(f, fp,  x0::Number; kwargs...)

    x   = copy(x0)
    tol = eps( eltype(real(x)) )
    D = [k => v for (k,v) in kwargs]
    xtol    = get(D, :xtol   , 100*tol)
    xtolrel = get(D, :xtolrel, tol    )
    ftol    = get(D, :ftol   , 100*tol)

    maxeval   = get(D, :maxeval  , 100)
    maxfneval = get(D, :maxfneval, 2000)

    F = ZeroFunction3(f, fp, f,
                      [x;],
                      f(x),
                      newton_update,
                      :initial, :na,
                      0, 1,
                      xtol, ftol, xtolrel)

    _findzero(F, maxeval; maxfneval=maxfneval)
end
"""

Implementation of Newton's method: `x_n1 = x_n - f(x_n)/ f'(x_n)`

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function=D(f)` -- derivative of `f`. Defaults to automatic derivative

* `x0::Real` -- initial guess

Keyword arguments:

* `ftol`. Stop iterating when |f(xn)| <= max(1, |xn|) * ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + max(1, |xn|) * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `maxfneval`. Stop iterating if more than this many function calls, throw error.

* `verbose::Bool=false` Set to `true` to see trace.

"""
# newton(f, fp, x0) now accepts a Complex guess
newton(f, fp, x0::Number; kwargs...) = fn_template(newtonmeth, f, fp, x0; kwargs...)
newton(f::Function, x0::Real; kwargs...) =  newton(f, D(f), float(x0); kwargs...)
newton(p::Poly, x0::Real; kwargs...) = newton(convert(Function, p), float(x0); kwargs...)

# Halley's method (cubic convergence)
function halley_update(F)
    xn = F.x[end]
    fxn, fpxn, fppxn = F.fxn, F.fp(xn), F.fpp(xn)

    xn1 = xn - 2fxn*fpxn / (2*fpxn*fpxn - fxn * fppxn)

    F.fxn = F.f(xn1)
    F.fncalls += 3
    push!(F.x, xn1)
end

function halleymeth(f, fp, fpp, x0::Real, args...;
                 kwargs...)

    x = float(x0)
    D = [k => v for (k,v) in kwargs]    
    xtol    = get(D, :xtol, 100*eps(eltype(x)))
    xtolrel = get(D, :xtolrel, eps(eltype(x)))
    ftol    = get(D, :ftol, 100*eps(eltype(x)))
    
    maxeval = get(D, :maxeval, 100)
    maxfneval = get(D, :maxfneval, 2000)
    
    F = ZeroFunction3(f, fp, fpp,
                      [x;],
                      f(x),
                      halley_update,
                      :initial, :na,
                      0, 1,
                      xtol, ftol, xtolrel)

    _findzero(F, maxeval; maxfneval=maxfneval)
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
halley(f, fp, fpp, x0::Real; kwargs...) = fn_template(halleymeth, f, fp, fpp, x0; kwargs...)
halley(f::Function, x0::Real; kwargs...) = halley(f, D(f), D2(f), float(x0); kwargs...)
halley(p::Poly, x0::Real; kwargs...) = halley(convert(Function, p), float(x0); kwargs...)
