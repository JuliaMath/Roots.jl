##################################################
## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley
## Historic, we have derivative free versions of similar order


## Newton
type Newton <: UnivariateZeroMethod
end

function init_state{T}(method::Newton, fs, x0::T)
    state = UnivariateZeroState(x0,
                                convert(T, typemax(real(x0))), # newton allows complex values, so a bit fussy
                                fs.f(x0),
                                fs.f(x0),
                                Nullable{Vector{T}}(),
                                0,
                                1,
                                false,
                                false,
                                false,
                                false,
                                "")
    state
end

function update_state{T}(method::Newton, fs, o::UnivariateZeroState{T})
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

function newton_setup{T}(f, fp, x0::T; 
                         xabstol=eps(), xreltol=eps(),
                         abstol=eps(), reltol=eps(),
                         maxevals=100, maxfnevals=100,
                         verbose::Bool=false,
                         kwargs...)

    prob = UnivariateZeroProblem(Newton(), f, fp, x0,  Nullable{Vector{T}}())
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol, maxevals, maxfnevals, verbose) 

    prob, options
end

find_zero(f, x0, method::Newton; kwargs...) = newton(f, x0; kwargs...)

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
function newton(f, fp, x0; kwargs...)
    prob, options = newton_setup(f, fp, float(x0); kwargs...)
    find_zero(prob,  Newton(), options)
end

newton(f, x0; kwargs...) = newton(f, D(f), x0; kwargs...)


## Halley


immutable SecondDerivative
    f
    fp
    fpp
end


type Halley <: UnivariateZeroMethod
end

function update_state{T}(method::Halley, fs, o::UnivariateZeroState{T})
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

function halley_setup{T <: AbstractFloat}(f, fp, fpp, x0::T;
                                          xabstol=eps(T), xreltol=eps(T),
                                          abstol=eps(T), reltol=eps(T),
                                          maxevals=100, maxfnevals=100, verbose::Bool=false)


    prob = UnivariateZeroProblem(Newton(), f, fp, fpp, x0,  Nullable{Vector{T}}())
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol, maxevals, maxfnevals, verbose) 

    prob, options
    
end

find_zero(f, x0, method::Halley; kwargs...) = halley(f, x0; kwargs...)


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
function halley(f, fp, fpp, x0; kwargs...)
    prob, options = halley_setup(f, fp, fpp, float(x0); kwargs...)
    find_zero(prob, Halley(), options)
end


halley(f, x0; kwargs...) = halley(f, D(f), D(f,2), float(x0); kwargs...)

