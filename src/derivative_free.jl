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

## A type to hold a function that we want to find the zero of.
abstract ZeroFunction 
type ZeroFunction1{T<:AbstractFloat} <: ZeroFunction
    f
    x::Vector{T}
    fxn::Vector{T}
    update::Function
    state::Symbol
    how::Symbol
    iterations::Int
    fncalls::Int
    ftol::T
    xtol::T
    xtolrel::T
    β::Float64
end
    
function Base.writemime(io::IO, ::MIME"text/plain", F::ZeroFunction)
    if F.state == :converged
        print(io, "xn=$(F.x[end]), iterations=$(F.iterations), fncalls=$(F.fncalls)")
    else
        print(io, "convergence failed, iterations=$(F.iterations), fncalls=$(F.fncalls)")
    end
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
function verbose_output(out::ZeroFunction)
    println("Steps = $(out.iterations), function calls = $(out.fncalls); steps:")
    xs = copy(out.x)
    x0 = shift!(xs)
    println(x0)
    while length(xs) > 0
        x1 = shift!(xs)
        printdiff(x1, x0)
        x0 = x1
    end
    println("")
end


## Algorithm to find a zero
function _findzero(f::ZeroFunction,
                   maxeval::Int= 100;
                   maxfneval::Int=2000,
                   kwargs...
                   )

    ## trivial check
    if abs(f.fxn[end]) <= max(1, abs(f.x[end])) * f.ftol
        f.state = :converged
        f.how = :ftol
        return(f)
    end
    
    f.state = :running
    
    while f.iterations < maxeval

        f.iterations += 1
        f.update(f)
        f.state != :running && return(f)

        xn = f.x[end]
        fxn = f.fxn[end]

        xtol = f.xtol +  max(1, abs(xn)) * f.xtolrel
        ftol = max(1, abs(xn)) * f.ftol
        
        if length(f.x) >= 2
            xn_1  = f.x[end-1]
            if abs(xn - xn_1) <= xtol
                ## possible convergence, but still expect fxn to be small
                if abs(fxn) > sqrt(ftol)
                    f.state = :failed
                    f.how = :xtol_not_ftol
                else
                    f.state = :converged
                    f.how = :xtol
                end
                break
            end
        end
        
        if abs(fxn) <= ftol
            f.state = :converged
            f.how = :ftol
            break
        end
        if f.fncalls >= maxfneval
            f.state = :failed
            f.how = :maxfneval
            break
        end
    end
    
    if f.state == :running
        if f.iterations >= maxeval
            f.state = :failed
            f.how = :maxeval
        
        else
            f.state = :huh
        end
    end

    f
end

## Basic structure of the functions (orders 16, 8, 5, 2, 1)
function _function_template_(updatefn, f, x0; kwargs...)
    x = map(float,[x0;])
    
    D = [k => v for (k,v) in kwargs]
    xtol    = get(D, :xtol,     100*eps(eltype(x)))
    xtolrel = get(D, :xtolrel,  eps(eltype(x)))
    ftol    = get(D, :ftol,     100*eps(eltype(x)))
    maxeval = get(D, :maxeval,  100)
    beta    = get(D, :beta,     1.0)

    F = ZeroFunction1(f,
                      x,
                      map(f,x),
                      updatefn,
                      :na, :na,
                      0, length(x),
                      ftol, xtol, xtolrel,
                      float(beta))
    
    
    _findzero(F, maxeval; kwargs...)
end


## check if abs(fx) <= tol
function check_ftol(fx, x, F)
    ftol = max(1,abs(x)) * F.ftol
    
    if abs(fx) > ftol
        false
    else
        push!(F.x, x)
        F.state = :converged
        F.how = :ftol
        true
    end
end

## check secant approximation
function check_secant(s, x, F)
    if abs(s) > max(1,abs(x)) * F.ftol
        false
    else
        F.state = :failed
        F.how = :zero_division
        true
    end
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
function thukral16_update(F)
    xn = F.x[end]
    fxn = F.fxn[end]

    check_ftol(fxn, xn, F) && return(F)

    wn = xn + 1/F.β * fxn  ## steffensen step
    fwn = F.f(wn)
    F.fncalls += 1
    check_ftol(fwn, wn, F) && return(F)

    secxnwn = (fxn - fwn) / (xn - wn)
    check_secant(secxnwn, xn, F) && return(F)
    
    yn = xn - fxn / secxnwn
    fyn = F.f(yn)
    F.fncalls += 1
    check_ftol(fyn, yn, F) && return(F)

    secxnyn = (fxn - fyn) / (xn - yn)
    secwnyn = (fwn - fyn) / (wn - yn)
    
    
    u3, u4 = fyn/fxn, fyn/fwn
    phi1, phi2, phi3 = 1/(1 - u4), 1 + u4 + u4^2, secxnwn / secwnyn

    check_secant(secxnyn, yn, F) && return(F)
    zn = yn - phi3 * fyn / secxnyn
    fzn = F.f(zn)
    F.fncalls += 1
    check_ftol(fzn, zn, F) && return(F)

    u1, u2 = fzn/fxn, fzn/fwn
    eta = 1/(1 + 2u3*u4^2)/(1 - u2)

    secynzn = (fyn - fzn) / (yn - zn)
    secxnzn = (fxn - fzn) / (xn - zn)
    appsec = secynzn - secxnyn + secxnzn

    check_secant(appsec, zn, F) && return(F)
    an = zn - eta * fzn / appsec
    fan = F.f(an)
    F.fncalls += 1
    check_ftol(fan, an, F) && return(F)

    u5, u6 = fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 + u2^2*u3 + 3u1*u4^2*(u3^2 - u4^2) / secxnyn
    secynan_1 = (yn - an)/(fyn - fan)
    secznan_1 = (zn - an)/(fzn -fan)
    xn1 = zn  - sigma * secynzn * fan * secynan_1 * secznan_1 

    F.fxn[end] = F.f(xn1)
    F.fncalls += 1
    push!(F.x, xn1)
end




## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
function thukral8_update(F)

    xn = F.x[end]
    fxn = F.fxn[end]
    
    if abs(fxn) <= F.ftol
        push!(F.x, xn)
        return(F)
    end

    wn  = xn + 1/F.β * fxn  # steffensen step
    fwn = F.f(wn)
    F.fncalls += 1
    check_ftol(fwn, wn, F) && return(F)

    secwnxn = (fwn - fxn) / fxn
    check_secant(secwnxn, xn, F) && return(F)
    
    yn = xn - fxn /secwnxn
    fyn = F.f(yn)
    F.fncalls += 1    
    check_ftol(fyn, yn, F) && return(F)

    
    ϕ = 1/(1 - fyn/fwn)
    secxnyn = (fxn - fyn) / (xn - yn)
    check_secant(secxnyn, yn, F) && return(F)
    
    zn = yn - ϕ * fyn / secxnyn
    fzn = F.f(zn)
    F.fncalls += 1
    check_ftol(fzn, zn, F) && return(F)

    secynzn = (fzn - fyn) / (zn - yn)
    secxnzn = (fzn - fxn) / (zn-xn)
    secxnynzn = secynzn - secxnyn + secxnzn
    check_secant(secxnynzn, zn, F) && return(F)
    
    ω = 1/(1 - fzn/fwn)
    ξ = 1 - 2 * fyn * fyn * fyn / (fwn * fwn * fxn)
    xn1 = zn - ω * ξ * fzn / secxnynzn

    F.fxn[end] = F.f(xn1)
    F.fncalls += 1
    push!(F.x, xn1)
end

## http://www.naturalspublishing.com/files/published/ahb21733nf19a5.pdf
## A New Fifth Order Derivative Free Newton-Type Method for Solving Nonlinear Equations
## Manoj Kumar, Akhilesh Kumar Singh, and Akanksha Srivastava
## Appl. Math. Inf. Sci. 9, No. 3, 1507-1513 (2015)
function kss5_update(F)

    xn = F.x[end]
    fxn = F.fxn[end]
    check_ftol(fxn, xn, F)
    
    wn = xn + 1/F.β * fxn # steffensen step
    fwn = F.f(wn)
    F.fncalls +=1
    check_ftol(fwn, wn, F) && return(F)

    secwnxn = (fwn - fxn) / fxn
    check_secant(secwnxn, xn, F) && return(F)
    
    yn = xn - fxn / secwnxn
    fyn = F.f(yn)
    F.fncalls +=1
    check_ftol(fyn, yn, F) && return(F)

    zn = xn - (fxn + fyn) / secwnxn
    fzn = F.f(zn)
    F.fncalls +=1    
    check_ftol(fzn, zn, F) && return(F)

    secxnyn = (fxn - fyn) / (xn - yn)
    secwnyn = (fwn - fyn) / (wn - yn)
    check_secant(secxnyn, zn, F) && return(F)
    check_secant(secwnyn, zn, F) && return(F)
    
    xn1 = zn - fzn * secwnxn / secxnyn / secwnyn

    F.fxn[end] = F.f(xn1)
    F.fncalls += 1
    push!(F.x, xn1)
end

## Order 2 methods
## http://en.wikipedia.org/wiki/Steffensen's_method
function steffensen_update(F)
    xn = F.x[end]
    fxn = F.fxn[end]

    Δf = F.f(xn+ 1/F.β * fxn) - fxn
    F.fncalls += 1

    secxn = Δf/fxn
    check_secant(secxn, xn, F) && return(F)
    xn1 = xn - fxn / secxn

    F.fxn[end] = F.f(xn1)
    F.fncalls += 1
    push!(F.x, xn1)
end

## order 1
function secmeth_update(F)
    a,b = F.x[end-1:end]
    fb,fa = F.fxn[end], F.fxn[end-1]
    
    x = b - (b-a) * fb/(fb - fa)

    F.fxn[(end-1):end] =  [fb, F.f(x)]
    F.fncalls += 1
    push!(F.x, x)
end



thukral16(f,x0; kwargs...)    = _function_template_(thukral16_update,  f, x0; kwargs...)
thukral8(f, x0; kwargs...)    = _function_template_(thukral8_update,   f, x0; kwargs...)
kss5(f, x0; kwargs...)        = _function_template_(kss5_update,       f, x0; kwargs...)
steffensen(f, x0; kwargs...)  = _function_template_(steffensen_update, f, x0; kwargs...)
secmeth(f, x0, x1; kwargs...) = _function_template_(secmeth_update,    f, [float(x0), float(x1)]; kwargs...)



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
        out = thukral16(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, kwargs...)
    elseif order == 8
        out = thukral8(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, kwargs...)        
    elseif order == 5
        out = kss5(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, kwargs...)                
    elseif order == 2
        out = steffensen(f, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, kwargs...)
    elseif order == 1
        x1 = x0 + 1e-4
        out = secmeth(f, x1, x0; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, kwargs...)
    else
      throw(ArgumentError())
    end

    if verbose
        verbose_output(out)
    end
        
    if out.state == :converged
        out.x[end]
    else
        ## failed
        steps = out.iterations > 1 ? "steps" : "step"
        steps_msg = "$(out.iterations) $steps"
        msg = ""
        if out.how == :zero_division
            msg = " Attempted division by 0."
        elseif out.how == :xtol_not_ftol
            msg = " Tolerance in Δx is small, but f(xn) is not within tolerance."
        elseif out.how == :maxfneval
            steps_msg = "$(out.fncalls) function calls"
        end
        throw(ConvergenceFailed("Failed to converge in $steps_msg.$msg Last value was $(out.x[end])."))
    end
end

