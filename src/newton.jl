##################################################
## Classical derivative-based, iterative, root-finding algorithms: Newton and Halley
## Historic, we have derivative free versions of similar order


## Newton
type Newton <: UnivariateZeroMethod
end

function callable_function(method::Newton, f::Tuple)
    length(f) == 1 && return FirstDerivative(f[1], D(f[1]))
    FirstDerivative(f[1], f[2])
end
callable_function(method::Newton, f::Any) = FirstDerivative(f, D(f))

function init_state{T}(method::Newton, fs, x0::T)
    state = UnivariateZeroState(x0,
                                x0 + typemax(Int),
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

function update_state{T}(method::Newton, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)
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
function derivative_free_setup{T<:AbstractFloat}(method::Newton, fs::CallableFunction, x0::Union{T, Complex{T}};
                                  bracket=Nullable{Vector{T}}(),
                                  xabstol=zero(T), xreltol=zero(T),
                                  abstol=4*eps(T), reltol=4*eps(T),
                                  maxevals=40, maxfnevals=typemax(Int),
                                  verbose::Bool=false)
    x = float(x0)

    if isa(x, Complex)
        bracket = Nullable{Vector{Complex{T}}}()     # bracket makes no sense for complex input, but one is expected
    elseif !isa(bracket, Nullable)
        bracket = Nullable(convert(Vector{T}, bracket))
    end

    prob = UnivariateZeroProblem(fs, x, bracket)
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol,  maxevals, maxfnevals, verbose)
    prob, options
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
newton(f, x0; kwargs...) = find_zero(f, x0, Newton(); kwargs...)
newton(f, fp, x0; kwargs...) = find_zero((f, fp), x0, Newton(); kwargs...)


## Halley



type Halley <: UnivariateZeroMethod
end

function callable_function(method::Halley, f::Tuple)
    length(f) == 1 && return SecondDerivative(f[1], D(f[1]), D(f[1],2))
    length(f) == 2 && return SecondDerivative(f[1], f[2], D(f[2],1))
    SecondDerivative(f[1], f[2], f[3])
end
callable_function(method::Halley, f) = SecondDerivative(f, D(f), D(f, 2))


function update_state{T}(method::Halley, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)
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
