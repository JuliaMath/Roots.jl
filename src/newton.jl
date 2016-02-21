## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley

## Newton
function newton_itr(f, fp, x0; 
                    xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
                    maxsteps=100, maxfnevals=100)
    update =  (o) -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]
        fpxn = o.fp(xn)
        incfn(o)

        xn1 = xn - fxn / fpxn
        fxn1 = o.f(xn1)
        incfn(o)

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end

    x, fx = promote(float(x0), f(float(x0)))
    out = ZeroType(f, fp, nothing, update, [x], [fx],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)


    out
end



function newton_method{T}(f, fp, x0::T;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    o = newton_itr(f, fp, x0;
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)


    for x in o
        nothing
    end
    verbose && verbose_output(o)
    o.xn[end]
end
newton_method{T<:AbstractFloat}(f, x::T; kwargs...) = newton_method(f, D(f), x; kwargs...)

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
newton{T<:Number}(f, fp, x0::T; kwargs...) = newton_method(f, fp, float(x0); kwargs...)
newton{T<:Real}(f, x0::T; kwargs...) = newton(f, D(f), float(x0); kwargs...)



## Halley's method (cubic convergence)
function halley_itr(f, fp, fpp, x0::Real;
                   xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
                   maxsteps=100, maxfnevals=100)

    update = (o) -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]
        fpxn = o.fp(xn)
        fppxn = o.fpp(xn)

        xn1 = xn - 2fxn*fpxn / (2*fpxn*fpxn - fxn * fppxn)
        fxn1 = o.f(xn1)

        incfn(o, 3)
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end

    x, fx = promote(float(x0), f(float(x0)))
    out = ZeroType(f, fp, fpp, update, [x], [fx],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end

function halley_method{T<:AbstractFloat}(f, fp, fpp, x0::T;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    o = halley_itr(f, fp, fpp, float(x0);
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)

    for x in o
        nothing
    end
    verbose && verbose_output(o)
    o.xn[end]
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
halley(f,fp, fpp, x::Number; kwargs...) = halley_method(f, fp, fpp, float(x); kwargs...)
halley(f, fp, x::Real; kwargs...) = halley(f, fp, D(fp), x; kwargs...)
halley(f, x::Real; kwargs...) = halley(f, D(f), D(f,2), x; kwargs...)
