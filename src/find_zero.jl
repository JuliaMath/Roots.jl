## Framework for setting up an iterative problem for finding a zero
## TODO
## * a graphic of trace when verbose=true?


# 
# In McNamee & Pan (DOI:10.1016/j.camwa.2011.11.015 there are a number of
# results on efficiencies of a solution, (1/d) log_10(q)
# Those implemented here are:
# quadratic cut (Muller) .265 (in a42)
# Newton() newton = .1505   or 1/2log(2)
# Order1() secant method .20 (1/1 * log(1.6)
# FalsePostion(12) Anderson-Bjork [.226, .233]
# FalsePostion(3) (King?) .264
# A42() 0.191 but convergence guaranteed
# Order8() 8th order 4 steps: .225 (log10(8)/4
# Order16() 16th order 5 steps .240
# Order5(): 5th order, 4 steps. 0.1747



# A zero is found by specifying:
# the method to use <: AbstractUnivariateZeroMethod
# the function(s) <: CallableFunction
# the initial state through a value for x either x, [a,b], or (a,b) <: AbstractUnivariateZeroState
# the options (e.g., tolerances) <: UnivariateZeroOptions

# The minimal amount needed to add a method, is to define a Method and an update_state method.

### Methods
abstract type AbstractUnivariateZeroMethod end
abstract type AbstractBisection <: AbstractUnivariateZeroMethod end
abstract type AbstractSecant <: AbstractUnivariateZeroMethod end


### States    
abstract type  AbstractUnivariateZeroState end
mutable struct UnivariateZeroState{T,S} <: AbstractUnivariateZeroState where {T,S}
    xn1::T
    xn0::Union{Missing, T}
    fxn1::S
    fxn0::Union{Missing, S}
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


# initialize state for most methods
function init_state(method::Any, fs, x)

    x1 = float(x)
    fx1 = fs(x1); fnevals = 1

    
    state = UnivariateZeroState(x1, missing, 
                                fx1, missing, 
                                0, fnevals,
                                false, false, false, false,
                                "")
    state
end



### Options
mutable struct UnivariateZeroOptions{Q,R,S,T}
    xabstol::Q
    xreltol::R
    abstol::S
    reltol::T
    maxevals::Int
    maxfnevals::Int
    verbose::Bool
end

# Allow for override of default tolerances. Useful, say, for methods like bisection
function _map_tolerance_arguments(d, xatol, xrtol, atol, rtol)
    xatol = get(d, :xabstol, xatol)
    xrtol = get(d, :xreltol, xrtol)
    atol = get(d, :abstol, atol)
    atol = get(d, :reltol, rtol)
    xatol, xrtol, atol, rtol
end

function init_options(::Any,
                      state;
                      xatol=missing,
                      xrtol=missing,
                      atol=missing,
                      rtol=missing,
                      maxevals::Int=40,
                      maxfnevals::Int=typemax(Int),
                      verbose::Bool=false,
                      kwargs...)

    ## Where we set defaults
    x1 = real(oneunit(float(state.xn1)))
    fx1 = real(oneunit(float(state.fxn1)))

    ## map old tol names to new
    ## deprecate in future
    xatol, xrtol, atol, rtol = _map_tolerance_arguments(Dict(kwargs), xatol, xrtol, atol, rtol)
    
    # assign defaults when missing
    options = UnivariateZeroOptions(ismissing(xatol) ? zero(x1) : xatol, # units of x
                                    ismissing(xrtol) ?  eps(x1/oneunit(x1)) : xrtol,  # unitless
                                    ismissing(atol)  ?  4.0 * eps(fx1) : atol,  # units of f(x)
                                    ismissing(rtol)  ?  4.0 * eps(fx1/oneunit(fx1)) : rtol, # unitless
                                    maxevals, maxfnevals,
    verbose)    

    options
end

### Functions
abstract type CallableFunction end

## It is faster the first time a function is used if we do not
## parameterize this. (As this requires less compilation) It is slower
## the second time a function is used. This seems like the proper
## tradeoff.  If it a case where the same function is being used many
## times, this function would be helpful
##
## function _find_zero(f, x0, method::Roots.AbstractUnivariateZeroMethod;kwargs...)
##     state = Roots.init_state(method, f, float.(x0))
##     options = Roots.init_options(method, state; kwargs...)
##     find_zero(method, f, options, state)
##  end    
struct DerivativeFree <: CallableFunction 
    f
end

struct FirstDerivative <: CallableFunction
    f
    fp
end

struct SecondDerivative <: CallableFunction
    f
    fp
    fpp
end

(F::DerivativeFree)(x::Number) = F.f(x)
(F::FirstDerivative)(x::Number) = F.f(x)
(F::SecondDerivative)(x::Number) = F.f(x)

(F::DerivativeFree)(x::Number, n::Int)  = F(x, Val{n})
(F::FirstDerivative)(x::Number, n::Int)  = F(x, Val{n})
(F::SecondDerivative)(x::Number, n::Int)  = F(x, Val{n})

(F::DerivativeFree)(x::Number, ::Type{Val{1}}) = D(F.f)(x)
(F::FirstDerivative)(x::Number, ::Type{Val{1}}) = F.fp(x)
(F::SecondDerivative)(x::Number, ::Type{Val{1}}) = F.fp(x)

(F::DerivativeFree)(x::Number, ::Type{Val{2}}) = D(F.f, 2)(x)
(F::FirstDerivative)(x::Number, ::Type{Val{2}}) = D(F.f, 2)(x)
(F::SecondDerivative)(x::Number, ::Type{Val{2}}) = F.fpp(x)


function callable_function(@nospecialize(fs))
    if isa(fs, Tuple)
        length(fs)==1 && return DerivativeFree(fs[1])
        length(fs)==2 && return FirstDerivative(fs[1],fs[2])
        return SecondDerivative(fs[1],fs[2],fs[3])
    end
    DerivativeFree(fs)
end
    


   
    
## has UnivariateZeroProblem converged?
## allow missing values in isapprox
function _isapprox(a, b, rtol, atol, lambda=missing, relaxed=false)
    _isapprox(Val{ismissing(a) || ismissing(b)}, float(a), float(b), rtol, atol,float(lambda),relaxed)
end
### missing data so not approx
_isapprox(::Type{Val{true}}, a, b, rtol, atol,lambda, relaxed) = false

function _isapprox(::Type{Val{false}}, a, b, rtol, atol, lambda, relaxed)

    if !ismissing(lambda)
        rtol *= max(one(lambda), abs(lambda/oneunit(lambda))) 
    end

    if relaxed
        tol = cbrt(max(atol/oneunit(atol), rtol))
        abs(a - b)/oneunit(a) <= tol
    else
        isapprox(a, b, rtol=rtol, atol=atol)
    end
end

function assess_convergence(method::Any, state, options)

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
        state.message = "NaN produced by algorithm."
        return true
    end
    
    if isinf(fxn1)
        state.convergence_failed = true
        state.message = "Inf produced by algorithm."
        return true
    end

    # f(xstar) ≈ xstar * f'(xstar)*eps(), so we pass in lambda
    if  _isapprox(fxn1, zero(fxn1), options.reltol, options.abstol, abs(xn1))
        state.f_converged = true
        return true
    end

    if _isapprox(xn1, xn0,  options.xreltol, options.xabstol)
        # Heuristic check that f is small too in unitless way
        λ = max(one(real(xn1)), abs(xn1/oneunit(xn1)))
        if _isapprox(fxn1, zero(fxn1), options.reltol, options.abstol, λ, true)
            state.x_converged = true
            return true
        end
    end


    if state.stopped
        if state.message == ""
            error("no message? XXX debug this XXX")
        end
        return true
    end

    return false
end

function show_trace(state, xns, fxns, method)
    converged = state.x_converged || state.f_converged
    
    println("Results of univariate zero finding:\n")
    if converged
        println("* Converged to: $(xns[end])")
        println("* Algorithm: $(method)")
        println("* iterations: $(state.steps)")
        println("* function evaluations: $(state.fnevals)")
        state.x_converged && println("* stopped as x_n ≈ x_{n-1} using atol=xatol, rtol=xrtol")
        state.f_converged && state.message == "" && println("* stopped as |f(x_n)| ≤ max(δ, max(1,|x|)⋅ϵ) using δ = atol, ϵ = rtol")
        state.message != "" && println("* Note: $(state.message)")
    else
        println("* Convergence failed: $(state.message)")
        println("* Algorithm $(method)")
    end
    println("")
    println("Trace:")
    
    itr, offset =  0:(lastindex(xns)-1), 1
    for i in itr
        x_i,fx_i, xi, fxi = "x_$i", "f(x_$i)", xns[i+offset], fxns[i+offset]
        println(@sprintf("%s = % 18.16f,\t %s = % 18.16f", x_i, float(xi), fx_i, float(fxi)))
    end
    println("")
    
    
end


"""

    find_zero(fs, x0, method; kwargs...)

Interface to one of several methods for find zeros of a univariate function.



# Initial starting value

For most methods, `x0` is a scalar value indicating the initial value in the iterative procedure. Values must be a subtype of `Number` and have methods for `float`, `real`, and `oneunit` defined. 

May also be a bracketing interval, specified as a tuple or a vector. A bracketing interval, (a,b), is one where f(a) and f(b) have different signs.

# Specifying a method

A method is specified to indicate which algorithm to employ:

* There are methods for bisection where a bracket is specified: `Bisection`

* There are methods for guarded bisection where a bracket is specified: `FalsePosition`

* There are several derivative-free methods: cf. `Order0`, `Order1` (secant method), `Order2` (Steffensen), `Order5`, `Order8`, and `Order16`, where the number indicates the order of the convergence.

* There are some classical methods where a derivative is assumed or computed using `ForwardDiff`: `Newton`, `Halley`. (The are not exported, so they need qualification, e.g., `Roots.Newton()`.

For more detail, see the help page for each method (e.g., `?Order5`).

If no method is specified, the default method depends on `x0`:

* If `x0` is a scalar, the default is the slower, but more robust `Order0` method.

* If `x0` is a tuple or vector indicating a *bracketing* interval, the `Bisection` method is use. (this method specializes on floating point values, but otherwise uses an algorithm of Alefeld, Potra, and Shi.)

# Specifying the function 

The function(s) are passed as the first argument. 

For the few methods that use a derivative (`Newton`, `Halley`, and optionally `Order5`)
a tuple of functions is used. For methods requiring a derivative and
second derivative, a tuple of three functions is used. If the
derivative functions are not specified, automatic differentiation via
the `ForwardDiff` package will be employed (for `Newton` and `Halley`).

# Optional arguments (tolerances, limit evaluations, tracing)

* `xatol` - absolute tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `xrtol` - relative tolerance for `x` values. Passed to `isapprox(x_n, x_{n-1})`
* `atol`  - absolute tolerance for `f(x)` values. Passed to `isapprox(f(x_n), zero(f(x_n))`
* `rtol`  - relative tolerance for `f(x)` values. Passed to `isapprox(f(x_n), zero(f(x_n))`
* `maxevals`   - limit on maximum number of iterations 
* `maxfnevals` - limit on maximum number of function evaluations
* `verbose` - if  `true` a trace of the algorithm will be shown on successful completion.

# Convergence

For most methods there are several heuristics used for convergence:

* if f(x_n) ≈ 0, using the tolerances `atol` and `rtol`, convergence is declared

* if x_n ≈ x_{n-1}, using the tolerances `xatol` and `xrtol`, *and* f(x_n) ≈ 0 with a relaxed tolerance then convergence is declared.

* if the algorithm has an issue (say a value of NaN appears) *and* f(x_n) ≈ 0 with a relaxed tolerance then convergence is declared, otherwise a failure to converge is declared

* if the number of iterations exceeds `maxevals` or the number of function evaluations exceeds `maxfnevals` a failure to converge is declared

* if x_n is `NaN` or f(x_n) is infinite  a failure to converge is declared


In general, with floating point numbers, convergence must be
understood as not an absolute statement. Even if mathematically x is
an answer the floating point realization, say xstar, may have
f(xstar) - f(x) = f(xstar) ≈ f'(x) ⋅ eps(x), so tolerances must be
appreciated, and at times specified.

For the `Bisection` methods, convergence is guaranteed, so the tolerances are set to be 0 by default.


# Examples:

```
# default methods
find_zero(sin, 3)  # use Order0()
find_zero(sin, (3,4)) # use Bisection()

# specifying a method
find_zero(sin, 3.0, Order2())              # Use Steffensen method
find_zero(sin, big(3.0), Order16())        # rapid convergence
find_zero(sin, (3, 4), FalsePosition())    # fewer function calls than Bisection(), in this case
find_zero(sin, (3, 4), FalsePosition(8))   # 1 or 12 possible algorithms for false position
find_zero((sin,cos), 3.0, Roots.Newton())  # use Newton's method
find_zero(sin, 3.0, Roots.Newton())        # use Newton's method with automatic f'(x)
find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley())  # use Halley's method

# changing tolerances
fn, x0, xstar = (x -> (2x*cos(x) + x^2 - 3)^10/(x^2 + 1), 3.0,  2.9806452794385368)
find_zero(fn, x0, Order2()) - xstar        # 0.011550654688925466
find_zero(fn, x0, Order2(), atol=0.0, rtol=0.0) # error: x_n ≉ x_{n-1}; just f(x_n) ≈ 0
fn, x0, xstar = (x -> (sin(x)*cos(x) - x^3 + 1)^9,        1.0,  1.117078770687451)
find_zero(fn, x0, Order2())                # 1.1122461983100858
find_zero(fn, x0, Order2(), maxevals=10)   # Roots.ConvergenceFailed: 26 iterations needed

# tracing output
find_zero(x->sin(x), 3.0, Order2(), verbose=true)   # 3 iterations
find_zero(x->sin(x)^5, 3.0, Order2(), verbose=true) # 23 iterations


```
"""
function find_zero(fs, x0, method::AbstractUnivariateZeroMethod; kwargs...)

    x = float.(x0)

    F = callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)

    find_zero(method, F, options, state)
    
end

find_zero(f, x0::T; kwargs...)  where {T <: Number}= find_zero(f, x0, Order0(); kwargs...)
find_zero(f, x0::Vector; kwargs...) = find_zero(f, x0, Bisection(); kwargs...)
find_zero(f, x0::Tuple; kwargs...) = find_zero(f, x0, Bisection(); kwargs...)

# Main method
function find_zero(M::AbstractUnivariateZeroMethod,
                   F,
                   options::UnivariateZeroOptions,
                   state::UnivariateZeroState{T,S}
                   )  where {T<:Number, S<:Number}


    
    # in case verbose=true
    if options.verbose
        if isa(M, AbstractSecant)
            xns, fxns = T[state.xn0, state.xn1], S[state.fxn0, state.fxn1]
        else
            xns, fxns = T[state.xn1], S[state.fxn1]
        end
    end

    while true
        
        val = assess_convergence(M, state, options)

        if val
            if state.stopped && !(state.x_converged || state.f_converged)
                ## stopped is a heuristic, there was an issue with an approximate derivative
                ## say it converged if pretty close, else say convergence failed.
                ## (Is this a good idea?)
                xstar, fxstar = state.xn1, state.fxn1

                λ = max(one(real(xstar)), abs(xstar/oneunit(xstar)))
                if _isapprox(fxstar, zero(fxstar), options.reltol, options.abstol, λ, true)
                    msg = "Algorithm stopped early, but |f(xn)| < ϵ^(1/3), where ϵ depends on xn, rtol, and atol"
                    state.message = state.message == "" ? msg : state.message * "\n\t" * msg
                    state.f_converged = true
                else
                    state.convergence_failed = true
                end
            end
                
            if state.x_converged || state.f_converged
                options.verbose && show_trace(state, xns, fxns, M)
                return state.xn1
            end

            if state.convergence_failed
                options.verbose && show_trace(state, xns, fxns, M)
                throw(ConvergenceFailed("Stopped at: xn = $(state.xn1)"))
                return state.xn1
            end
        end

        update_state(M, F, state, options)

        if options.verbose
            push!(xns, state.xn1)
            push!(fxns, state.fxn1)
        end

    end
end

