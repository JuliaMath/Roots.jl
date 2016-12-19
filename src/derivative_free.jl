## TOD0
## add no-op maxevals, maxfnevals, ... to bisection methods

# container for callable objects
immutable DerivativeFree
    f
end

immutable FirstDerivative
    f
    fp
end

## method names
abstract UnivariateZeroMethod

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

# generic initialization. Modify as necessary for a method
# we use typemax as a sentinel, as we check |xn-xn_1| and |f(xn) - f(xn_1)| for smallness
# using Nullable would be better, but currently this is a hassle
function init_state{T}(method::Any, fs, x0::T, bracket)
    fx0 = fs.f(x0)
    S = eltype(fx0)
    state = UnivariateZeroState(x0,
                                x0 + typemax(Int),
                                fx0,
                                fx0,
                                isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, bracket)),
                                0,
                                1,
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

## basic container
type UnivariateZeroProblem
    fs
    x0
    bracket::Nullable
end

## override for diferent methods if need be.
function UnivariateZeroProblem{T}(method::UnivariateZeroMethod, f, x0, bracket::Nullable{Vector{T}})
    UnivariateZeroProblem(DerivativeFree(f), x0, bracket)
end
function UnivariateZeroProblem{T}(method::UnivariateZeroMethod, f, fp, x0::T, bracket::Nullable{Vector{T}})
    UnivariateZeroProblem(FirstDerivative(f, fp), x0, bracket)
end
function UnivariateZeroProblem{T}(method::UnivariateZeroMethod, f, fp, fpp, x0::T, bracket::Nullable{Vector{T}})
    UnivariateZeroProblem(SecondDerivative(f,fp, fpp), x0, bracket)
end



## set up a derivative free problem
## output is run through find_zero
function derivative_free_setup{T<: AbstractFloat}(method, f, x0::Union{T,Vector{T}};
                                                  bracket=Nullable{Vector{T}}(),
                                                  xabstol=zero(T), xreltol=zero(T),
                                                  abstol=4*eps(T), reltol=4*eps(T),
                                                  maxevals=40, maxfnevals=typemax(Int),
                                                  verbose::Bool=false)

    prob = UnivariateZeroProblem(method, f, x0, bracket)
    options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol,  maxevals, maxfnevals, verbose)
    prob, options
end

## has UnivariateZeroProblem converged?
function assess_convergence(fs, state, options)

    xn0, xn1 = state.xn0, state.xn1
    fxn0, fxn1 = state.fxn0, state.fxn1

    
    if (state.x_converged || state.f_converged)
        return true
    end

    if state.steps > options.maxevals
        state.convergence_failed = true
        state.message = "too many steps taken"
        return true
    end

    if state.fnevals > options.maxfnevals
        state.convergence_failed = true
        state.message = "too many function evaluations taken"
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
    
    λ = norm(xn1, Inf) 
    
    if norm(fxn1) <= max(options.abstol, λ * options.reltol)
        state.f_converged = true
        return true
    end

    if norm(xn1 - xn0) <= max(options.xabstol, λ * options.xreltol)
        state.x_converged = true
        return true
    end


    if state.stopped
        if state.message == ""

            error("no message?")
        end
        return true
    end

    return false
end

function show_trace(state, xns, fxns, method)
    converged = state.x_converged || state.f_converged
    
    println("Results of univariate zero finding:\n")
    if converged
        println("* Converged to $(xns[end])")
        println("* Algorithm $(method)")
        println("* iterations: $(state.steps)")
        println("* function evaluations: $(state.fnevals)")
        state.x_converged && println("* stopped as |x_n - x_{n-1}| < ϵ")
        state.f_converged && println("* stopped as |f(x_n)| < ϵ")
    else
        println("Convergence failed: $(state.message)")
    end
    println("")
    println("Trace:")
    for i in eachindex(xns)-1
        x_i,fx_i, xi, fxi = "x_$i", "f(x_$i)", xns[i+1], fxns[i+1]
        println(@sprintf("%s = %18.15f,\t %s = %18.15f", x_i, float(xi), fx_i, float(fxi)))
    end
    println("")
end



function find_zero(prob::UnivariateZeroProblem,  method::UnivariateZeroMethod, options::UnivariateZeroOptions)
    fs = prob.fs
    state = init_state(method, fs, prob.x0, prob.bracket)
    find_zero(method, fs, state, options)
end

function find_zero(method::UnivariateZeroMethod, fs, state::UnivariateZeroState, options::UnivariateZeroOptions)
    xns, fxns = [state.xn1], [state.fxn1]

    while true
        if !isnull(state.bracket)
            m,M = get(state.bracket)
            if (state.xn1 < m || state.xn1 > M) || state.stopped
                # do bisection step, update bracket
                c = m + (M-m)/2
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

        val = assess_convergence(fs, state, options)

        if options.verbose
            push!(xns, state.xn1)
            push!(fxns, state.fxn1)
        end


        if val
            if state.stopped
                ## stopped is a heuristic, there was an issue with an approximate derivative
                ## say it converged if pretty close, else say convergence failed.
                xstar, fxstar = state.xn1, state.fxn1
                if abs(fxstar) <= (options.abstol)^(2/3)
                    state.f_converged = true
                else
                    state.convergence_failed = true
                end
            end
                
            if state.x_converged || state.f_converged
                options.verbose && show_trace(state, xns, fxns, method)
                return state.xn1
            end

            if state.convergence_failed
                options.verbose && show_trace(state, xns, fxns, method)
                options.verbose && warn(state.message)
                throw(ConvergenceFailed(""))
            end
        end

        update_state(method, fs, state)
    end
end

"""

Find a zero of a univariate function using one of several different methods.

Positional arugments:

* `f` a function or callable object. Most methods are derivative free. Some (`Newton`, `Halley`) may have derivative(s) computed using the `ForwardDiff` pacakge.

* `x0` an initial starting value. Typically a scalar, but may be an array for bisection methods. The value `float(x0)` is passed on.

* `method` one of several methods, see below.

Keyword arguments:

* `xabstol=zero()`: declare convergence if |x_n - x_{n-1}| <= max(xabstol, norm(x_n) * xreltol)

* `xreltol=eps()`:

* `abstol=zero()`: declare convergence if |f(x_n)}| <= max(abstol, norm(x_n) * reltol)

* `reltol`:

* `bracket`: Optional. A bracketing interval for the sought after
  root. If given, a hybrid algorithm may be used where bisection is
  utilized for steps that would go out of bounds.

* `maxevals::Int=40`: stop trying after `maxevals` steps

* `maxfnevals::Int=typemax(Int)`: stop trying after `maxfnevals` function evaluations

* `verbose::Bool=false`: If `true` show information about algorithm

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
function find_zero(f, x0, method::Any; kwargs...)
    T = method
    prob, options = derivative_free_setup(T, f, map(float,x0); kwargs...)
    find_zero(prob, T, options)
end



## old interface
function derivative_free{T <: AbstractFloat}(f, x0::T;
                                             ftol::Real = 10*eps(T),
                                             xtol::Real =  zero(T),
                                             xreltol::Real = 4.0 * eps(T),
                                             maxevals::Int = 30,
                                             verbose::Bool=false,
    order::Int=0,  # 0, 1, 2, 5, 8 or 16
    kwargs...      # maxfnevals,
                         )
    
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

    find_zero(f, x0, method;
              xabstol=xtol, xreltol=xreltol,
              abstol = ftol, reltol = ftol,
              maxevals=maxevals, maxfnevals=4*maxevals,
              verbose=verbose)

end



##################################################
## Helpers for the various methods
## issue with approx derivative
isissue(x) = (x == 0.0) || isnan(x) || isinf(x)


"""
heuristic to get a decent first step with Steffensen steps
"""
function steff_step{T}(x::T, fx)
    thresh =  max(1, norm(x)) * sqrt(eps(T)) # max(1, sqrt(abs(x/fx))) * 1e-6
    norm(fx) <= thresh ? fx : sign(fx) * thresh
end

function guarded_secant_step(alpha, beta, falpha, fbeta)
    fp = (fbeta - falpha) /  (beta - alpha) 
    Δ = fbeta / fp

    if norm(Δ) >= 100 * norm(alpha - beta) # guard runaway
        Δ = sign(Δ) * 100 * norm(alpha - beta)
    end
    
    beta - Δ, isissue(Δ)
    
end


## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    num, den = fb-fa, b - a
    num==0 && den == 0 && return Inf, true
    out = num / den
    out, isissue(out)    
end

## use f[y,z] - f[x,y] + f[x,z] to approximate
function _fbracket_diff(a,b,c, fa, fb, fc)
    x1, state = _fbracket(b, c, fb,  fc)
    x2, state = _fbracket(a, b, fa,  fb)
    x3, state = _fbracket(a, c, fa,  fc)
    out = x1 - x2 + x3
    out, isissue(out)
end


## use f[a,b] * f[a,c] / f[b,c]
function _fbracket_ratio(a, b, c, fa, fb, fc)
    x1,_ = _fbracket(b, c, fb, fc)
    x2,_ = _fbracket(a, b, fa, fb)
    x3,_ = _fbracket(a, c, fa, fc)
    out = (x2 * x3) / x1
    out, isissue(out)    
end

## Bisection
type Bisection end
function find_zero(f, x0, method::Bisection; xabstol=eps(typeof(float(x0[1]))), verbose=false, kwargs...)
    find_zero_bisection(f, float(x0[1]), float(x0[2]); xtol=xabstol, verbose=verbose) 
end



## Order0 and Secant are related
type Order0 <: UnivariateZeroMethod
end
type Secant <: UnivariateZeroMethod
end
typealias Order1 Secant

    

function init_state{T}(method::Union{Order0, Order1}, fs, x0::T, x1::T, bracket)

    fx0, fx1 = fs.f(x0), fs.f(x1)
    S = eltype(fx0)
    
    state = UnivariateZeroState(
                                x0, x1,
                                fx0, fx1,
                                isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, bracket)),
                                0,
                                2,
                                false,
                                false,
                                false,
                                false,
                                "")
    state
end
function init_state{T}(method::Union{Order0, Order1}, fs, x::T, bracket)
    # guess for x1
    if isa(x, Vector)
        x0, x1 = x[1:2]
    else
        x0 = float(x)
        fx0 = fs.f(x0)
        stepsize = max(1/100, min(abs(fx0), abs(x0/100)))
        x1 = x0 + stepsize
    end
    init_state(method, fs, x0, x1, bracket)
end


## order 0
function update_state{T}(method::Order0, fs, o::UnivariateZeroState{T})
    f = fs.f
    alpha, beta = o.xn0, o.xn1
    falpha, fbeta = o.fxn0, o.fxn1
    S = eltype(falpha)
    

    if sign(falpha) * sign(fbeta) < 0.0
        # use bisection
        o.xn1 = find_zero(f, [float(alpha), float(beta)], Bisection())
        incsteps(o, 53)
        incfn(o, 54)
        o.fxn1 = f(o.xn1)
        o.message = "Used bisection"        
        o.x_converged = true        
        return nothing
    end

    gamma, issue = guarded_secant_step(alpha, beta, falpha, fbeta)
    if issue
        o.message = "error with guarded secant step"
        o.stopped = true
        return nothing
    end

    fgamma = f(gamma)
    if sign(fgamma)*sign(fbeta) < 0
        m,M = min(gamma, beta), max(gamma,beta)
        o.xn1 = find_zero(f, map(float, [m, M]), Bisection())
        ## or pass in state...
        o.fxn1 = f(o.xn1)
        incsteps(o, 53)
        incfn(o, 54)
        o.message = "Used bisection"
        o.x_converged = true
        return nothing
    end

    if norm(fgamma) <= norm(fbeta)
        o.xn0, o.xn1 = beta, gamma
        o.fxn0, o.fxn1 = fbeta, f(gamma)
        return nothing
    end
    
    ctr = 0
    while true
        ## quadratic step
        ctr += 1
        # quadratic_step. Put new gamma at vertex of parabola through alpha, beta, (old) gamma
        num = falpha * (beta^2 - gamma^2) + fbeta * (gamma^2 - alpha^2) + fgamma * (alpha^2 - beta^2)
        den = falpha*(beta - gamma) + fbeta*(gamma - alpha) + fgamma*(alpha - beta)
        gamma = num / (2 * den)
        
        if ctr >= 5
            o.stopped = true
            o.message = "Failed to improve with quadratic step"
            return nothing
        end
        fgamma = f(gamma)
        incfn(o)
        if norm(fgamma) < norm(fbeta)
            o.xn0, o.xn1 = beta, gamma
            o.fxn0, o.fxn1 = fbeta, fgamma
            return nothing
        end
    end

    # failed to improve
    o.stopped = true
    o.message("failure to improve")
    return nothing
end



    
## Secant
function update_state{T}(method::Secant, fs, o::UnivariateZeroState{T})

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    S = eltype(fxn0)

    fp, issue = _fbracket(xn0, xn1, fxn0, fxn1)
    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end
    
    xn2::T = xn1 -  fxn1 / fp
    fxn2::S = fs.f(xn2)
    incfn(o)

    o.xn0, o.xn1 = xn1, xn2
    o.fxn0, o.fxn1 = fxn1, fxn2

    incsteps(o)

    nothing
end



# function secant_method{T <: AbstractFloat}(f, x0::T, x1::T;  bracket=Nullable{Vector{T}}(),
#                                            xabstol=zero(T), xreltol=zero(T),
#                                            abstol=4*eps(T), reltol=4*eps(T),
#                                            maxevals=40, maxfnevals=typemax(Int),
#                                            verbose::Bool=false)
#     fs = DerivativeFree(f)
#     prob = UnivariateZeroProblem(Secant(), f, [x0,x1], bracket)
#     options = UnivariateZeroOptions(xabstol, xreltol, abstol, reltol, maxevals, maxfnevals, verbose)


#     find_zero(prob, Secant(), options)
# end
# secant_method(f, x0::Real, x1::Real; kwargs...) = secant_method(f, float(x0), float(x1); kwargs...)
secant_method(f, x0::Real, x1::Real; kwargs...) = find_zero(f, map(float, [x0,x1]), Order1(); kwargs...)




### Steffensen
type Steffensen <: UnivariateZeroMethod
end
typealias Order2 Steffensen

function update_state{T}(method::Steffensen, fs, o::UnivariateZeroState{T})

    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn) 
    incfn(o)
        
    fp, issue = _fbracket(xn, wn, fxn, fwn)
        
    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end
        
    xn1::T = xn - fxn / fp
        
        
    fxn1::S = fs.f(xn1)
    incfn(o)
        
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1


    nothing
end

steffenson(f, x0; kwargs...) = find_zero(f, x0, Steffensen(); kwargs...)

##################################################
        


## A New Fifth Order Derivative Free Newton-Type Method for Solving Nonlinear Equations
## Manoj Kumar, Akhilesh Kumar Singh, and Akanksha Srivastava
## Appl. Math. Inf. Sci. 9, No. 3, 1507-1513 (2015)
## http://www.naturalspublishing.com/files/published/ahb21733nf19a5.pdf
type Order5 <: UnivariateZeroMethod
end

function update_state{T}(method::Order5, fs::DerivativeFree, o::UnivariateZeroState{T})
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)
    
    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
    incfn(o)
    
    fp, issue = _fbracket(xn, wn, fxn, fwn)
    if issue
        o.xn0, o.xn1 = xn, wn
        o.fxn0, o.fxn1 = fxn, fwn
        o.stopped  = true
        return
    end
    
    yn::T = xn - fxn / fp
    fyn::S = fs.f(yn)
    incfn(o)
    

    zn::T = xn - (fxn + fyn) / fp
    fzn::S = fs.f(zn)
    incfn(o)

    fp1 = _fbracket(xn, yn, fxn, fyn)[1] * _fbracket(wn, yn, fwn, fyn)[1] / fp

    if isissue(fp1)
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped = true
        return
    end
    
    xn1::T = zn  - fzn  / fp1
    fxn1::S = fs.f(xn1)
    incfn(o)
        
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end

## If we have a derivative
function update_state{T}(method::Order5, fs::FirstDerivative, o::UnivariateZeroState{T})
    incsteps(o)

    xn, fxn = o.xn1, o.fxn1
    S = eltype(fxn)

    fpxn = fs.fp(xn)
    incfn(o)

    if isissue(fpxn)
        o.stopped  = true
        return
    end
    
    yn = xn - fxn / fpxn
    fyn, fpyn = fs.f(yn), fs.fp(yn)
    incfn(o, 2)

    if isissue(fpyn)
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped  = true
        return
    end

    
    zn = xn  - (fxn + fyn) / fpxn
    fzn = fs.f(zn)
    incfn(o, 1)

    xn1 = zn - fzn / fpyn
    fxn1 = fs.f(xn1)
    incfn(o, 1)
    
    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end

kss5(f, x0; kwargs...) = kss5(f, D(f), x0; kwargs...)
function kss5(f, fp, x0; kwargs...)
    fs = FirstDerivative(f, fp)
    T = typeof(float(x0))
    state = init_state(Order5(), fs, x0, Nullable{Vector{T}}())
    options = UnivariateZeroOptions(eps(T), eps(T), eps(T), eps(T), 40, 200, true)
    zp = UnivariateZeroProblem(fs, state, options)

    find_zero(zp, Order5())
end


##################################################

## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
type Order8 <: UnivariateZeroMethod
end


function update_state{T}(method::Order8, fs, o::UnivariateZeroState{T})
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)
    
    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
    incfn(o)

    fp, issue = _fbracket(xn, wn, fxn, fwn)
    issue && return (xn, true)

    if issue
        o.stopped = true
        return
    end
    
    yn::T = xn - fxn / fp
    fyn::S = fs.f(yn)
    incfn(o)

    fp, issue = _fbracket(yn, xn, fyn, fxn)
    if issue
        o.xn0,o.xn1 = xn, yn
        o.fxn0,o.fxn1 = fxn, fyn
        o.stopped = true
        return
    end

    phi = (1 + fyn / fwn)           # pick one of options
    zn =  yn - phi * fyn / fp
    fzn::S = fs.f(zn)
    incfn(o)
        
    fp, issue =  _fbracket_diff(xn, yn, zn, fxn, fyn, fzn) 
    if issue
        o.xn0,o.xn1 = xn, zn
        o.fxn0,o.fxn1 = fxn, fzn
        o.stopped = true
        return
    end

    w::T = 1 / (1 - fzn/fwn)
    xi::T = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))
    
    xn1::T = zn - w * xi * fzn / fp
    fxn1::S = fs.f(xn1)
    incfn(o)

    o.xn0,o.xn1 = xn, xn1
    o.fxn0,o.fxn1 = fxn, fxn1
            
    nothing
end



##################################################

## 16th order, derivative free root finding algorithm
## http://article.sapub.org/10.5923.j.ajcam.20120203.08.html
## American Journal of Computational and Applied Mathematics
## p-ISSN: 2165-8935    e-ISSN: 2165-8943
## 2012;  2(3): 112-118
## doi: 10.5923/j.ajcam.20120203.08
## New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations
## R. Thukral
## Research Centre, 39 Deanswood Hill, Leeds, West Yorkshire, LS17 5JS, England
## from p 114 (17)
type Order16 <: UnivariateZeroMethod
end

function update_state{T}(method::Order16, fs, o::UnivariateZeroState{T})
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)
    
    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
    incfn(o)

    fp, issue = _fbracket(xn, wn, fxn, fwn)
    if issue
        o.xn0, o.xn1 = xn, wn
        o.fxn0, o.fxn1 = fxn, fwn
        o.stopped = true
        return
    end
        
    yn::T = xn - fxn / fp
    fyn::S = fs.f(yn)
    incfn(o)
    

    fp, issue = _fbracket(xn, yn, fxn, fyn)
    phi = _fbracket(xn, wn, fxn, fwn)[1] / _fbracket(yn, wn, fyn, fwn)[1]
    if issue
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped = true
        return
    end

    
    zn::T = yn - phi * fyn / fp
    fzn::S = fs.f(zn)
    incfn(o)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn
    eta = 1 / (1 + 2*u3*u4^2) / (1 - u2)
    if issue
        o.xn0, o.xn1 = xn, zn
        o.fxn0, o.fxn1 = fxn, fzn
        o.stopped = true
        o.message = "Approximate derivative failed"
        return
    end        

    an = zn - eta * fzn / fp
    fan = fs.f(an)
    incfn(o)

    fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
    sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
        u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/_fbracket(xn,yn, fxn, fyn)[1]
    
    if issue
        o.xn0, o.xn1 = xn, an
        o.fxn0, o.fxn1 = fxn, fan
        o.stopped = true
        o.message("Approximate derivative failed")        
        return
    end  
        
    xn1 = an - sigma * fan / fp
    fxn1 = fs.f(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


