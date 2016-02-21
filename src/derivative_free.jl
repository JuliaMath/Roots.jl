## Various derivative-free iterative algorithms to find a zero

## Use iterator for solving.
## allow for solving along the lines of:
## ```
## for xn in itr
##    println(xn)
## end
## ```
type ZeroType{T, S}
    f
    fp
    fpp
    
    update

    xn::Vector{T}
    fxn::Vector{S}

    xtol
    xtolrel
    ftol
    state
    
    cnt
    maxcnt
    fevals
    maxfevals
end
incfn(o::ZeroType,i::Int=1) = (o.fevals = o.fevals + i)
inccnt(o::ZeroType, i::Int=1) = (o.cnt = o.cnt + i)


function Base.start{T}(o::ZeroType{T})
    (o.xn[end], o.fxn[end])
end

function Base.next{T,S}(o::ZeroType{T,S}, state)
    o.update(o)
    o.cnt = o.cnt + 1
    (o.xn[end], o.fxn[end]), (o.xn[end], o.fxn[end])
end

## done function centralizes the stopping rules
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
    ftol = lambda * o.ftol
    xtol = o.xtol +  lambda * o.xtolrel

    norm(o.fxn[end]) < ftol && return true
    if length(o.xn) > 1
        if (norm(o.xn[end] - o.xn[end-1])) <= xtol && (norm(o.fxn[end]) < sqrt(ftol))
            return true
        end
    end
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

## issue with approx derivative
isissue(x) = (x == 0.0) || isnan(x) || isinf(x)


"""
heuristic to get a decent first step with Steffensen steps
"""
function steff_step(x, fx)
    thresh =  max(1, sqrt(abs(x/fx))) * 1e-6
    abs(fx) <= thresh ? fx : sign(fx) * thresh
end


## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    out = (fb - fa) / (b - a)
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



##################################################
## Iterators

# iterator for secant function
function secant_itr(f, x0::Real, x1::Real; xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(), maxsteps=100, maxfnevals=100)
    update = (o) -> begin
        xn_1, xn = o.xn[(end-1):end]
        fxn_1, fxn = o.fxn[(end-1):end]

        fp, issue = _fbracket(xn, xn_1, fxn, fxn_1)

        
        xn1 = xn - fxn  / fp
        fxn1 = o.f(xn1)
        incfn(o)
        
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end
    x, fx = promote(float(x1), f(float(x1)))
    out = ZeroType(f, nothing, nothing, update,
                   [float(x0), x], [f(float(x0)), fx],
                   xtol, xtolrel, ftol, :not_converged,
                   0, maxsteps, 2, maxfnevals)

    out
end


function steffensen_itr{T<:AbstractFloat}(f, x0::T;
                   xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
                   maxsteps=100, maxfnevals=100)
    update = (o) -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]
        
        wn = xn + steff_step(xn, fxn)
        fwn = f(wn)
        incfn(o)
        
        fp, issue = _fbracket(xn, wn, fxn, fwn)
        
        if issue
            o.state = :converged
            return
        end
        
        xn1 = xn - fxn / fp
        
        
        fxn1 = o.f(xn1)
        incfn(o)
        
        push!(o.xn, xn1)
        push!(o.fxn, fxn1)

    end

    x = float(x0)
    x, fx = promote(x, f(x))
    out = ZeroType(f, nothing, nothing, update, [x], [fx],
                   xtol, xtolrel, ftol, :not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end


## http://www.naturalspublishing.com/files/published/ahb21733nf19a5.pdf
## A New Fifth Order Derivative Free Newton-Type Method for Solving Nonlinear Equations
## Manoj Kumar, Akhilesh Kumar Singh, and Akanksha Srivastava
## Appl. Math. Inf. Sci. 9, No. 3, 1507-1513 (2015)
function kss5_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + steff_step(xn, fxn)
        fwn = o.f(wn)
        incfn(o)

        fp, issue = _fbracket(xn, wn, fxn, fwn)
        if issue
            o.state = :converged
            return
        end
        
        yn = xn - fxn / fp
        fyn = o.f(yn)
        incfn(o)
        
        zn = xn - (fxn + fyn) / fp ## not a step in thukral algorithms
        fzn = o.f(zn)
        incfn(o)
        fp, issue = _fbracket_ratio(yn, xn, wn, fyn, fxn, fwn)

        if issue || isinf(fzn)
            push!(o.xn, yn)
            push!(o.fxn, fyn)
            o.state = :converged
            return
        end
        
        xn1 = zn  - fzn  / fp
        fxn1 = o.f(xn1); incfn(o)
        

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    x, fx = promote(float(x0), f(float(x0)))
    out = ZeroType(f, nothing, nothing, update, [x], [fx],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
function thukral8_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + steff_step(xn, fxn)
        fwn = o.f(wn)
        incfn(o)

        fp, issue = _fbracket(xn, wn, fxn, fwn)
        issue && return (xn, true)

        if issue
            o.state = :converged
            return
        end

        yn = xn - fxn / fp
        fyn = o.f(yn)
        incfn(o)

        fp, issue = _fbracket(yn, xn, fyn, fxn)
        if issue
            push!(o.xn, yn); push!(o.fxn, fyn)
            o.state = :converged
            return
         end

        phi = (1 + fyn / fwn)           # pick one of options
        zn =  yn - phi * fyn / fp
        fzn = o.f(zn)
        incfn(o)
        
        fp, issue =  _fbracket_diff(xn, yn, zn, fxn, fyn, fzn) 
        if issue
            push!(o.xn, zn)
            push!(o.fxn, fzn)
            o.state = :converged
            return
         end

        w = 1 / (1 - fzn/fwn)
        xi = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))

        xn1 = zn - w * xi * fzn / fp
        fxn1 = o.f(xn1)
        incfn(o)

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    x, fx = promote(float(x0), f(float(x0)))    
    out = ZeroType(f, nothing, nothing, update, [x], [fx],
                   xtol, xtolrel, ftol,:not_converged, 
                   0, maxsteps, 1, maxfnevals)

    out
end


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
function thukral16_itr(f, x0::Real;
             xtol=4*eps(), xtolrel=4*eps(), ftol=4*eps(),
             maxsteps=100, maxfnevals=100)

    update = o -> begin
        xn = o.xn[end]
        fxn = o.fxn[end]

        wn = xn + steff_step(xn, fxn)
        fwn = o.f(wn)
        incfn(o)

        fp, issue = _fbracket(xn, wn, fxn, fwn)
        if issue 
            o.state = :converged
            return
         end
        
        yn = xn - fxn / fp
        fyn = o.f(yn)
        incfn(o)
        

        fp, issue = _fbracket(xn, yn, fxn, fyn)
        phi = _fbracket(xn, wn, fxn, fwn)[1] / _fbracket(yn, wn, fyn, fwn)[1]
        if issue
            push!(o.xn, yn); push!(o.fxn, fyn)
            o.state = :converged
            return
         end

        
        zn = yn - phi * fyn / fp
        fzn = o.f(zn)
        incfn(o)

        fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
        u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn
        eta = 1 / (1 + 2*u3*u4^2) / (1 - u2)
        if issue
            push!(o.xn, zn); push!(o.fxn, fzn)
            o.state = :converged
            return
        end        

        an = zn - eta * fzn / fp
        fan = o.f(an)
        incfn(o)

        fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
        u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
        sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
                 u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/_fbracket(xn,yn, fxn, fyn)[1]
        
        if issue
            push!(o.xn, an); push!(o.fxn, fan)
            o.state = :converged
            return
        end  
        
        xn1 = an - sigma * fan / fp
        fxn1 = o.f(xn1)
        incfn(o)

        push!(o.xn, xn1)
        push!(o.fxn, fxn1)
    end

    x, fx = promote(float(x0), f(float(x0)))    
    out = ZeroType(f, nothing, nothing, update, [x], [fx],
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

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol + abs(1, |xn|) * xtolrel. Checks that f(xn) is reasonably close.

* `xtolrel`. Stop iterating when |xn+1 - xn| <= xtol + abs(1, |xn|) * xtolrel. Checks that f(xn) is reasonably close.

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `maxfneval`. Stop iterating if more than this many function calls, throw error.

* `order`. One of 0, 1, 2, 5, 8, or 16. Specifies which algorithm to
use.

- Default is 0 for the slower, more robust, SOLVE function. 
- order 1 is a secant method
- order 2 is a Steffensen method
- order 5 uses a method of Kumar, Kumar Singh, and Srivastava
- order 8 From Thukral. Seems a bit more robust than the secant method and Steffensen method
- order 16 A higher-order method due to Thurkal. It may be faster when used with `Big` values.

* `verbose`. If true, will print number of iterations, function calls, and each step taken

* `kwargs...` passed on.

The `SOLVE` method has different stopping criteria.

The file test/test_fzero2 generates some timed comparisons

Because these are derivative free, they can be used with functions
defined by automatic differentiation
e.g., find critical point of f(x) = x^2
```
fzero(D(x -> x^2), 1)
```
"""

function derivative_free{T <: AbstractFloat}(f, x0::T;
                         ftol::Real = 10.0 * eps(x0),
                         xtol::Real =  4.0 * eps(x0),
                         xtolrel::Real = eps(x0),
                         maxeval::Int = 30,
                         verbose::Bool=false,
                         order::Int=0,  # 0, 1, 2, 5, 8 or 16
                         kwargs...      # maxfnevals,
                         )
    
    order == 0 && return(SOLVE(f, x0; ftol=ftol, maxeval=maxeval, verbose=verbose, kwargs...))
    _derivative_free(f, x0; order=order, ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval,
                     verbose=verbose, kwargs...)

end
function _derivative_free{T <: AbstractFloat}(f, x0::T;
                         ftol::Real = 10.0 * eps(x0),
                         xtol::Real =  4.0 * eps(x0),
                         xtolrel::Real = eps(x0),
                         maxeval::Int = 30,
                         verbose::Bool=false,
                         order::Int=0,  # 0, 1, 2, 5, 8 or 16
                         kwargs...      # maxfnevals, possible beta to control steffensen step
                         )
    if order == 16
        o = thukral16_itr(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval, kwargs...)
    elseif order == 8
        o = thukral8_itr(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval, kwargs...)        
    elseif order == 5
        o = kss5_itr(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)                
    elseif order == 2
        o = steffensen_itr(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)
    elseif order == 1
        x1 = x0 + steff_step(x0, f(x0))
        o = secant_itr(f, float(x1), x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxsteps=maxeval,  kwargs...)
    else
        throw(ArgumentError())
    end

    if done(o, start(o))
        verbose && println("Done before we started...")
        return(x0)::T
    end

    val = x0
    try
        for (x,fx) in o
            nothing
        end
        verbose && verbose_output(o)        
        o.xn[end]
    catch err
        verbose && verbose_output(o)                
        rethrow(err)
    end


end


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
function secant_method(f, x0::Real, x1::Real;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)

    x_0, x_1 = float(x0), float(x1)
    o = secant_itr(f, x_0, x_1;
                   xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                   maxsteps=maxsteps, maxfnevals=maxfnevals)

    for x in o
        nothing
    end

    verbose && verbose_output(o)        
    o.xn[end]
end

function steffensen_method{T<:AbstractFloat}(f, x0::T;
                       xtol=4*eps(), xtolrel=4*eps(), ftol=4eps(),
                       maxsteps::Int=100, maxfnevals=100,
                       verbose::Bool=false)
    o = steffensen_itr(f, x0; 
                  xtol=xtol, xtolrel=xtolrel, ftol=ftol,
                  maxsteps=maxsteps, maxfnevals=maxfnevals)

    for x in o
        nothing
    end
    verbose && verbose_output(o)
    o.xn[end]
end
steffensen(f, x0::Number, args...; kwargs...) = steffensen_method(f, float(x0), args...; kwargs...)

