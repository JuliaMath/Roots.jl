## Framework for setting up an iterative problem for finding a zero
## Different methods can implement:
## - A subtype of UnivariateZeroMethod (reqd)
## - UnivariateZeroProblem
## - callable_function
## - init_state
## - assess_convergence
## - udpate_state (reqd)
## TODO
## * a graphic of trace when verbose=true?

## method names are subtypes
@compat abstract type UnivariateZeroMethod end

# container for callable objects; not really necessary, but has some value.
@compat abstract type CallableFunction end
immutable DerivativeFree <: CallableFunction
    f
end

immutable FirstDerivative <: CallableFunction
    f
    fp
end


immutable SecondDerivative <: CallableFunction
    f
    fp
    fpp
end

## allows override for automatic derivatives, see Newton
function callable_function(m::UnivariateZeroMethod, f)
    !isa(f, Tuple) && return DerivativeFree(f)
    length(f) == 1 && return DerivativeFree(f[1])
    length(f) == 2 && return FirstDerivative(f[1], f[2])
    SecondDerivative(f[1], f[2], f[3])
end


## object to hold state
type UnivariateZeroState{T,S} 
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    bracket::Nullable{Vector{T}}
    steps::Int
    fnevals::Int
    stopped::Bool             # stopped, butmay not have converged
    x_converged::Bool         # converged via |x_n - x_{n-1}| < ϵ
    f_converged::Bool         # converged via |f(x_n)| < ϵ
    convergence_failed::Bool
    message::AbstractString
end

incfn(o::UnivariateZeroState, k=1)    = o.fnevals += k
incsteps(o::UnivariateZeroState, k=1) = o.steps += k

## generic initialization. Modify as necessary for a method, such as secant which uses both xn1 and xn0
## we use x0 + typemax(Int) as a sentinel. This could be a Nullable, but that
## is a bit more hassle
function init_state{T}(method::Any, fs, x0::T, bracket)
    fx0 = fs.f(x0); fnevals = 1

    if isa(bracket, Nullable)
        if !isnull(bracket)
            a,b = get(bracket)
            sign(fs.f(a)) * sign(fs.f(b)) > 0 && (warn(bracketing_error); throw(ArgumentError))
            fnevals += 2
        end
    else
        a,b = bracket[1,2]
        sign(fs.f(a)) * sign(fs.f(b)) > 0 &&  (warn(bracketing_error); throw(ArgumentError))
        fnevals += 2
    end
    
    S = eltype(fx0)
    state = UnivariateZeroState(x0,
                                x0 + typemax(Int),
                                fx0,
                                fx0,
                                isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, bracket)),
                                0,
                                fnevals,
                                false,
                                false,
                                false,
                                false,
                                "")
    state
end


## Options for convergence, reporting
type UnivariateZeroOptions{T}
    xabstol::T
    xreltol::T
    abstol::T
    reltol::T
    maxevals::Int
    maxfnevals::Int
    verbose::Bool
end


function univariate_zero_options{T}(args...;
                                    xabstol::T=zero(T),
                                    xreltol::T=zero(T),
                                    abstol::T=4*eps(T),
                                    reltol::T=4*eps(T),
                                    maxevals::Int=40,
                                    maxfnevals::Int=typemax(Int),
                                    verbose::Bool=false,
                                    kwargs...)

    ## adjust for old argument names
    kw = Dict(kwargs)
    z = zero(T)
    UnivariateZeroOptions(max(z, get(kw, :xtol, xabstol)),
                          max(z, get(kw, :xtolrel, xreltol)),
                          max(z, abstol),
                          max(z, get(kw, :ftol, reltol)),
                          max(0, get(kw, :maxeval, get(kw, :maxsteps, maxevals))),
                          max(0, get(kw, :maxfneval, maxfnevals)),
                          verbose)
end
        

## basic container
type UnivariateZeroProblem{T<:AbstractFloat}
    fs::CallableFunction
    x0::Union{T, Vector{T}, Complex{T}}
    bracket::Nullable
end

## frame the problem and the options
function derivative_free_setup(method::Any, fs::CallableFunction, x0; kwargs...)
    _derivative_free_setup(method, fs, x0; kwargs...)
end
function _derivative_free_setup{T<:AbstractFloat}(method::Any, fs::CallableFunction, x0::Union{T, Vector{T}};
                                                 bracket=Nullable{Vector{T}}(),
                                                 xabstol=zero(T), xreltol=zero(T),
                                                 abstol=4*eps(T), reltol=4*eps(T),
                                                 maxevals=40, maxfnevals=typemax(Int),
                                                 verbose::Bool=false,
                                                 kwargs...
    )
    x = map(float,x0)
    prob = UnivariateZeroProblem(fs, x, isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, bracket)))
    options = univariate_zero_options(;xabstol=xabstol,
                                      xreltol=xreltol,
                                      abstol=abstol,
                                      reltol=reltol,
                                      maxevals=maxevals,
                                      maxfnevals=maxfnevals,
                                      verbose=verbose,
                                      kwargs...)
    prob, options
end


## has UnivariateZeroProblem converged?
function assess_convergence(method::Any, fs, state, options)

    xn0, xn1 = state.xn0, state.xn1
    fxn0, fxn1 = state.fxn0, state.fxn1

    
    if (state.x_converged || state.f_converged)
        return true
    end

    if state.steps > options.maxevals
        state.stopped = true
        state.message = "too many steps taken."
        return true
    end

    if state.fnevals > options.maxfnevals
        state.stopped = true
        state.message = "too many function evaluations taken."
        return true
    end

    if isnan(xn1)
        state.convergence_failed = true
        state.message = "NaN produced by algorithm"
        return true
    end
    
    if isinf(fxn1)
        state.convergence_failed = true
        state.message = "Inf produced by algorithm"
        return true
    end

    λ = max(1, norm(xn1))
    
    if  norm(fxn1) <= max(options.abstol, λ * options.reltol)
        state.f_converged = true
        return true
    end

    if isapprox(xn1, xn0, rtol = options.xreltol, atol=options.xabstol) && norm(fxn1) <= cbrt(max(options.abstol, λ * options.reltol))
        state.x_converged = true
        return true
    end


    if state.stopped
        if state.message == ""
            error("no message? XXX debug this XXX")
        end
        return true
    end

    return false
end

function show_trace(fs, state, xns, fxns, method)
    converged = state.x_converged || state.f_converged
    
    println("Results of univariate zero finding:\n")
    if converged
        println("* Converged to: $(xns[end])")
        println("* Algorithm $(method)")
        println("* iterations: $(state.steps)")
        println("* function evaluations: $(state.fnevals)")
        state.x_converged && println("* stopped as x_n ≈ x_{n-1} using atol=xabstol, rtol=xreltol")
        state.f_converged && state.message == "" && println("* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = abstol, ϵ = reltol")
        state.message != "" && println("* Note: $(state.message)")
    else
        println("* Convergence failed: $(state.message)")
        println("* Algorithm $(method)")
    end
    println("")
    println("Trace:")
    
    itr, offset =  0:(endof(xns)-1), 1
    for i in itr
        x_i,fx_i, xi, fxi = "x_$i", "f(x_$i)", xns[i+offset], fxns[i+offset]
        println(@sprintf("%s = % 18.16f,\t %s = % 18.16f", x_i, float(xi), fx_i, float(fxi)))
    end
    println("")

    
end

### find_zero method has potentially many interfaces

## problem, method, options. This may be a main entry point
function find_zero(prob::UnivariateZeroProblem,  method::UnivariateZeroMethod, options::UnivariateZeroOptions)
    fs = prob.fs
    state = init_state(method, fs, prob.x0, prob.bracket)
    find_zero(method, fs, state, options)
end

## method, state, options. Could be used to switch  between methods
function find_zero(method::UnivariateZeroMethod, fs, state::UnivariateZeroState, options::UnivariateZeroOptions)

    # in case verbose=true
    if isa(method, AbstractSecant)
        xns, fxns = [state.xn0, state.xn1], [state.fxn0, state.fxn1]
    else
        xns, fxns = [state.xn1], [state.fxn1]
    end
    
    while true
        if  !isnull(state.bracket)
            m,M = get(state.bracket)
            if (state.xn1 < m || state.xn1 > M) || state.stopped
                # do bisection step, update bracket
                c = m + (M-m)/2
                # don't land on xn0
                if c == state.xn0
                    c = m + (M-m)/4
                end
                fxn1 = fs.f(c)
                state.xn1, state.fxn1 = c, fxn1
                if sign(fxn1) * sign(fs.f(m)) < 0
                    state.bracket = Nullable([m,c])
                else
                    state.bracket = Nullable([c,M])
                end
                state.stopped && (state.stopped = false)
                incfn(state)
            end
        end
        
        val = assess_convergence(method, fs, state, options)



        if val
            if state.stopped
                ## stopped is a heuristic, there was an issue with an approximate derivative
                ## say it converged if pretty close, else say convergence failed.
                ## (Is this a good idea?)
                xstar, fxstar = state.xn1, state.fxn1
                if abs(fxstar) <= (options.abstol)^(2/3)
                    msg = "Algorithm stopped early, but |f(xn)| < ϵ^(2/3), where ϵ = abstol"
                    state.message = state.message == "" ? msg : state.message * "\n\t" * msg
                    state.f_converged = true
                else
                    state.convergence_failed = true
                end
            end
                
            if state.x_converged || state.f_converged
                options.verbose && show_trace(fs, state, xns, fxns, method)
                return state.xn1
            end

            if state.convergence_failed
                options.verbose && show_trace(fs, state, xns, fxns, method)
                throw(ConvergenceFailed("Stopped at: xn = $(state.xn1)"))
            end
        end

        update_state(method, fs, state, options)

        if options.verbose
            push!(xns, state.xn1)
            push!(fxns, state.fxn1)
        end

    end
end

"""

Find a zero of a univariate function using one of several different methods.

Positional arugments:

* `f` a function, callable object, or tuple of same. A tuple is used
  to pass in derivatives, as desired. Most methods are derivative
  free. Some (`Newton`, `Halley`) may have derivative(s) computed
  using the `ForwardDiff` pacakge.

* `x0` an initial starting value. Typically a scalar, but may be an
  array for bisection methods. The value `float(x0)` is passed on.

* `method` one of several methods, see below.

Keyword arguments:

* `xabstol=zero()`: declare convergence if |x_n - x_{n-1}| <= max(xabstol, max(1, |x_n|) * xreltol)

* `xreltol=eps()`:

* `abstol=zero()`: declare convergence if |f(x_n)| <= max(abstol, max(1, |x_n|) * reltol)

* `reltol`:

* `bracket`: Optional. A bracketing interval for the sought after
  root. If given, a hybrid algorithm may be used where bisection is
  utilized for steps that would go out of bounds.

* `maxevals::Int=40`: stop trying after `maxevals` steps

* `maxfnevals::Int=typemax(Int)`: stop trying after `maxfnevals` function evaluations

* `verbose::Bool=false`: If `true` show information about algorithm and a trace.

Returns: 

Returns `xn` if the algorithm converges. If the algorithm stops, returns `xn` if 
|f(xn)| ≤ ϵ^(2/3), where ϵ = reltol, otherwise a `ConvergenceFailed` error is thrown.

Exported methods: 

`Bisection()`;
`Order0()` (heuristic, slow more robust);
`Order1()` (also `Secant()`);
`Order2()` (also `Steffensen()`);
`Order5()` (KSS);
`Order8()` (Thukral);
`Order16()` (Thukral);

Not exported:

`Secant()`, use `Order1()`
`Steffensen()` use `Order2()`
`Newton()` (use `newton()` function)
`Halley()` (use `halley()` function)

The order 0 method is more robust to the initial starting point, but
can utilize many more function calls. The higher order methods may be
of use when greater precision is desired.`


Examples:

```
f(x) = x^5 - x - 1
find_zero(f, 1.0, Order5())
find_zero(f, 1.0, Steffensen()) # also Order2()
```
"""
find_zero{T<:Number}(f, x0::Union{T,Vector{T}}, method::UnivariateZeroMethod; kwargs...) =
    find_zero(method, callable_function(method, f), x0; kwargs...)

## some defaults for methods
find_zero{T <: Number}(f, x0::T; kwargs...) = find_zero(f, x0, Order0(); kwargs...)
find_zero{T <: Number}(f, x0::Vector{T}; kwargs...) = find_zero(f, x0, Bisection(); kwargs...)


function find_zero(method::UnivariateZeroMethod, fs::CallableFunction, x0; kwargs...)
    x = float(x0)
    prob, options = derivative_free_setup(method, fs, x; kwargs...)
    find_zero(prob, method, options)
end





## old interface for fzero
## old keyword arguments (see ?fzero) handled in univariate_zero_options
function derivative_free{T <: AbstractFloat}(f, x0::T; order::Int=0,
                                             kwargs...)
    
    if order == 0
        method = Order0()
    elseif order == 1
        method = Order1()
    elseif order == 2
        method = Order2()
    elseif order == 5
        method = Order5()
    elseif order == 8
        method = Order8()
    elseif order == 16
        method = Order16()
    else
        warn("Invalid order. Valid orders are 0, 1, 2, 5, 8, and 16")
        throw(ArgumentError())
    end

    find_zero(f, x0, method; kwargs...)
end

