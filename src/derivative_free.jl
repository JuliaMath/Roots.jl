## Various higher order, derivative-free root finding algorithms
##
## These are basically based on an update step which is 3 or 4
## Newton's method calls of the type xn1 = xn - w f(xn)/f'(xn) with
## different choices for the weights and different derivative-free
## approximations for f'(xn).
##
## The 2nd order Steffensen method is faster at times than the 5, 8 or
## 16th order ones. The 5th and 8th seem more robust to
## initial condition and take fewer function calls than the 16th.

## Use iterator for solving

type ZeroType{T <: AbstractFloat, S <: AbstractFloat}
    f
    fp
    fpp
    
    update

    xn::Vector{T}
    fxn::Vector{S}

    xtol
    xtolrel
    ftol
    state::Symbol
    
    cnt
    maxcnt
    fevals
    maxfevals
end

function Base.start(o::ZeroType)
    (o.xn[end], o.fxn[end])
end

function Base.next(o::ZeroType, state)
    o.update(o)
    o.cnt = o.cnt + 1
    o.xn[end], (o.xn[end], o.fxn[end])
end


function Base.done(o::ZeroType, state)
    if o.state == :converged
        ## check for near convergence
        (norm(o.fxn[end]) <= sqrt(o.ftol)) && return true
        throw(ConvergenceFailed("Algorithm stopped with xn=$(o.xn[end]), f(xn)=$(o.fxn[end])"))
    end
        
    o.cnt > o.maxcnt && throw(ConvergenceFailed("Too many steps taken"))
    o.fevals > o.maxfevals  && throw(ConvergenceFailed("Too many function calls taken"))
    isnan(o.fxn[end]) && throw(ConvergenceFailed("NaN produced by algorithm"))
    isinf(o.xn[end]) && throw(ConvergenceFailed("Algorithm escaped to oo"))
                               
    # return turn if f(xn) \approx 0 or xn+1 - xn \approx 0
    lambda = max(1, abs(o.xn[end]))
    xtol = o.xtol +  lambda * o.xtolrel
    ftol = lambda * o.ftol

    norm(o.fxn[end]) < ftol && return true
    (length(o.xn)>1 && norm(o.xn[end] - o.xn[end-1]) <= xtol) && return true
    false
end


## printing utility, borrowed idea of showing change from previous from @stevengj. Only works at REPL
function printdiff(x1, x0, color=:red)
    s0, s1 = string(x0), string(x1)
    flag = true
    for i in 1:length(s1)
        if i <= length(s0) && s0[i:i] == s1[i:i]
            print_with_color(:red, s1[i:i])
        else
            print(s1[i:end])
            break
        end
    end
    print("\n")
end

# show output
function verbose_output(out::ZeroType)
    println("Steps = $(out.cnt), function calls = $(out.fevals); steps:")
    xs = copy(out.xn)
    x0 = shift!(xs)
    println(x0)
    while length(xs) > 0
        x1 = shift!(xs)
        printdiff(x1, x0)
        x0 = x1
    end
    println("")
end


##################################################
isissue(x) = isnan(x) || isinf(x)


function secant_itr(f, x0::Real, x1::Real; xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(), maxsteps=100, maxfnevals=100)
    update = (o) -> begin
        xn_1, xn = o.xn[(end-1):end]
        fxn_1, fxn = o.fxn[(end-1):end]
        
        xn1 = xn - fxn * (xn - xn_1) / (fxn - fxn_1)
        fxn1 = o.f(xn1); o.fevals = o.fevals + 1
        
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end
    out = ZeroType(f, nothing, nothing, update,
                   [float(x0), float(x1)], [float(f(x0)), float(f(x1))],
                   xtol, xtolrel, ftol, :not_converged,
                   0, maxsteps, 2, maxfnevals)

    out
end

function secant_method(f, x0::Real, x1::Real;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)
    
    o = secant_itr(f, float(x0), float(x1);
                   xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                   maxsteps=maxsteps, maxfnevals=maxfnevals)
    try
        collect(o)
        verbose && verbose_output(o)        
    catch e
        verbose && verbose_output(o)        
        rethrow(e)
    end
    o.xn[end]
end



function steffensen_itr(f, x0::Real;
                   xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
                   maxsteps=100, maxfnevals=100)
    update = (o) -> begin
        xn = o.xn[end]
        h = fxn = o.fxn[end]
        
        xn1 = xn - fxn / (o.f(xn + fxn) - fxn) * fxn;    o.fevals = o.fevals + 1
        fxn1 = o.f(xn1);  o.fevals = o.fevals + 1
        
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end
        
    out = ZeroType(f, nothing, nothing, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol, :not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end
function steffensen_method(f, x0::Real;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    o = steffensen_itr(f, float(x0); 
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)
    collect(o)
    verbose && verbose_output(o)
    o.xn[end]
end
steffensen(args...; kwargs...) = steffensen_method(args...; kwargs...)

## Newton
function newton_itr(f, fp, x0; 
                    xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
                    maxsteps=100, maxfnevals=100)
    update =  (o) -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]
        fpxn = o.fp(xn);  o.fevals = o.fevals + 1

        xn1 = xn - fxn / fpxn
        fxn1 = o.f(xn1);  o.fevals = o.fevals + 1

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end
        
    out = ZeroType(f, fp, nothing, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)


    out
end
function newton_method{T<:Number}(f, fp, x0::T;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    o = newton_itr(f, fp, float(x0);
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)

    collect(o)
    verbose && verbose_output(o)
    o.xn[end]
end
newton_method(f, x::Real; kwargs...) = newton_method(f, D(f), x; kwargs...)
newton(args...; kwargs...) = newton_method(args...; kwargs...)

## Halley
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

        o.fevals = o.fevals + 3
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end
        
    out = ZeroType(f, fp, fpp, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end
function halley_method(f, fp, fpp, x0::Real;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    o = halley_itr(f, fp, fpp, float(x0);
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)

    collect(o)
    verbose && verbose_output(o)
    o.xn[end]
end
halley_method(f::Function, fp::Function, x::Number; kwargs...) = halley_method(f, fp, fp', x; kwargs...)
halley_method(f::Function, x::Number; kwargs...) = halley_method(f, D(f), D(f,2), x; kwargs...)
halley(args...; kwargs...) = halley_method(args...; kwargs...)

function kss5_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + fxn
        fwn = o.f(wn); o.fevals = o.fevals + 1

        fwnxn_i = (wn - xn) / (fwn - fxn)
        if isissue(fwnxn_i)
            o.state = :converged
            return
        end
        
        yn = xn - fxn * fwnxn_i
        fyn = o.f(yn); o.fevals = o.fevals + 1
        
        zn = xn - (fxn - fyn) * fwnxn_i
        fzn = o.f(zn); o.fevals = o.fevals + 1
        
        fxnyn_i =  (xn - yn) / (fxn - fyn)
        fwnyn_i = (wn - yn) / (fwn - fyn)
        fwnxn = (fwn - fxn) / (wn - xn)
        
        if (isissue(fwnxn_i) | isissue(fwnyn_i))
            push!(o.xn, yn)
            push!(o.fxn, fyn)
            o.state = :converged
            return
        end
        
        xn1 = zn  - fzn  * fwnxn * fxnyn_i * fwnyn_i
        fxn1 = o.f(xn1); o.fevals = o.fevals + 1
        

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    out = ZeroType(f, nothing, nothing, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end


function thukral8_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + fxn
        fwn = o.f(wn); o.fevals = o.fevals + 1
        
        yn = xn - fxn * fxn / (fwn - fxn)
        fyn = o.f(yn); o.fevals = o.fevals + 1
        
        fynxn_1 = (yn - xn) / (fyn - fxn)
        if isissue(fynxn_1)
            push!(o.xn, yn); push!(o.fxn, fyn)
            o.state = :converged
            return
         end

        phi = (1 + fyn / fwn)           # pick one of options
        zn =  yn - phi * fynxn_1 * fyn
        fzn = o.f(zn);  o.fevals = o.fevals + 1
        
        fznyn = (fzn - fyn) / (zn - yn)
        fynxn = (fyn - fxn) / (yn - xn)
        fznxn = (fzn - fxn) / (zn - xn)
        fp_1  = 1 / (fznyn - fynxn + fznxn)
         if isissue(fp_1)
            push!(o.xn, zn); push!(o.fxn, fzn)
            o.state = :converged
            return
         end

        w = 1 / (1 - fzn/fwn)
        xi = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))

        xn1 = zn - w * xi * fp_1 * fzn
        fxn1 = o.f(xn1); o.fevals = o.fevals + 1

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    out = ZeroType(f, nothing, nothing, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end



function thukral16_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + fxn
        fwn = o.f(wn); o.fevals = o.fevals + 1

        fwnxn_1 = (wn - xn) / (fwn - fxn)
        if isissue(fwnxn_1) 
            o.state = :converged
            return
         end
        
        yn = xn - fxn * fwnxn_1
        fyn = o.f(yn); o.fevals = o.fevals + 1
        

        fxnyn_1 = (xn - yn) / (fxn - fyn)
        if isissue(fwnxn_1)
            push!(o.xn, yn); push!(o.fxn, fyn)
            o.state = :converged
            return
         end

        zn = yn - fyn * fxnyn_1
        fzn = o.f(zn); o.fevals = o.fevals + 1

        fznyn = (fzn - fyn) / (zn - yn)
        fynxn = (fyn - fxn) / (yn - xn)
        fznxn = (fzn - fxn) / (zn - xn)
        fp_1 = 1 / (fznyn - fynxn + fznxn)
        if isissue(fp_1)
            push!(o.xn, zn); push!(o.fxn, fzn)
            o.state = :converged
            return
        end        

        an = zn - fzn * fp_1
        fan = o.f(an); o.fevals = o.fevals + 1

        fynxn = (fyn - fxn) / (yn - xn)
        fynan_1 = (yn - an) / (fyn - fan)
        fznan_1 = (zn - an) / (fzn - fan)
        if isissue(fynan_1) | isissue(fznan_1)
            push!(o.xn, an); push!(o.fxn, fan)
            o.state = :converged
            return
        end  
        
        xn1 = an - fan * fynxn * fynan_1 * fznan_1
        fxn1 = o.f(xn1); o.fevals = o.fevals + 1

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    out = ZeroType(f, nothing, nothing, update, [float(x0)], [float(f(x0))],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end



"""
Main interface for derivative free methods

* `f` a scalar function `f:R -> R` or callable object. Methods try to find a solution to `f(x) = 0`.

* `x0` initial guess for zero. Iterative methods need a reasonable
initial starting point.

* `ftol`. Stop iterating when |f(xn)| <= max(1, |xn|) * ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + abs(1, |xn|) * xtolrel

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + abs(1, |xn|) * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `maxfneval`. Stop iterating if more than this many function calls, throw error.

* `order`. One of 0, 1, 2, 5, 8, or 16. Specifies which algorithm to
use.

- Default is 0 for the slower, more robust, SOLVE function. 
- order 1 is a secant method
- order 2 is a Steffensen method
- order 5 uses a method of Kumar, Kumar Singh, and Srivastava
- order 8 of Thukral is a bit more robust than the secant method and Steffensen method
- order 16 is a higher-order method due to Thurkal. It may be faster when used with `Big` values.

* `verbose`. If true, will print number of iterations, function calls, and each step taken

* `kwargs...` passed on. For orders 2, 5, 8, and 16 a value for `beta`
  can control initial Steffensen step by `1/beta*f(xn)`

The `SOLVE` method has different stopping criteria.

The file test/test_fzero2 generates some timed comparisons

Because these are derivative free, they can be used with functions
defined by automatic differentiation
e.g., find critical point of f(x) = x^2
```
fzero(D(x -> x^2), 1)
```
"""
function derivative_free(f, x0::Real;
                         ftol::Real = 10.0 * eps(float(x0)),
                         xtol::Real =  4.0 * eps(float(x0)),
                         xtolrel::Real = eps(float(x0)),
                         maxeval::Int = 30,
                         verbose::Bool=false,
                         order::Int=0,  # 0, 1, 2, 5, 8 or 16
                         kwargs...      # maxfnevals, possible beta to control steffensen step
                         )

    order == 0 && return(SOLVE(f, x0; ftol=ftol, maxeval=maxeval, verbose=verbose, kwargs...))

    
    if order == 16
        o = thukral16_itr(f, float(x0); ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval, kwargs...)
    elseif order == 8
        o = thukral8_itr(f, float(x0); ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval, kwargs...)        
    elseif order == 5
        o = kss5_itr(f, float(x0); ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)                
    elseif order == 2
        o = steffensen_itr(f, float(x0); ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)
    elseif order == 1
        x1 = x0 + 1e-4
        o = secant_itr(f, float(x1), float(x0); ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)
    else
      throw(ArgumentError())
    end

    
    collect(o)
    verbose && verbose_output(o)
    o.xn[end]
end

