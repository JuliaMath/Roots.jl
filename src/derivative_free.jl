## Various higher order, derivative free root finding algorithms
##
## These are basically based on an update step which is 3 or 4
## Newton's method calls of the type xn1 = xn - w f(xn)/f'(xn) with
## different choices for the weights and different derivative-free
## approximations for f'(xn).
##
## The 2nd order Steffensen method is faster at times than the 5, 8 or
## 16th order ones, but 8th is the default as it seems more robust to
## initial condition, as we replace an initial Steffensen step with a
## secant line approximation when f(x) is too large.

## some helpers
function secant(fa::Real, fb::Real, a::Real, b::Real)
    (fa - fb)/(a-b)
end

function approx_deriv(f::Function, fx::Real, x::Real, ftol=0.0)
    abs(fx) <= ftol && throw(StateConverged(x))
    (f(x + fx) - fx)/fx
end

function secant2(f, fz, fx, z, x, ftol=0.0) 
    if abs(x-z) == 0.0
        abs(fx) <= ftol && throw(StateConverged(x))
        throw(DomainError())
    end
    (approx_deriv(f, fx, x) - secant(fz, fx, z, x))/(x-z)
end

function F(f::Function, fx, fy, fz, x, y, z, ftol=0.0)
    secant(fz, fy, z, y)  + secant2(f, fz, fx, z, x, ftol) * (z-y)
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
function thukral_update16(f::Function, x0::Real, ftol::Real; kwargs...)

    xn = x0
    fxn = f(xn)
    abs(fxn) <= ftol && throw(StateConverged(xn))

    ## first step is of size fxn. If that is big, this whole thing has
    ## problems with convergence. Here we replace with a step based on
    ## the central difference
    if abs(fxn) > 1e-1
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end
    
    wn = xn + fxn
    fwn = f(wn)
    abs(fwn) <= ftol && throw(StateConverged(wn))


    secxnwn = secant(fxn, fwn, xn, wn)
    ((secxnwn == 0.0) | isnan(secxnwn)) && throw(ConvergenceFailed("Division by 0"))
    
    yn = xn - fxn/secxnwn
    fyn = f(yn)
    abs(fyn) <= ftol && throw(StateConverged(wn))

    secxnyn = secant(fxn, fyn, xn, yn)
    ((secxnyn == 0.0) | isnan(secxnyn)) && throw(ConvergenceFailed("Division by 0"))
    
    secwnyn = secant(fwn, fyn, wn, yn)
    ((secwnyn == 0.0) | isnan(secwnyn)) && throw(ConvergenceFailed("Division by 0"))
      
    u3, u4 = fyn/fxn, fyn/fwn
    phi1, phi2, phi3 = 1/(1 - u4), 1 + u4 + u4^2, secxnwn / secwnyn


    zn = yn - phi3 * fyn / secxnyn
    fzn = f(zn)
    abs(fzn) <= ftol && throw(StateConverged(zn))
    
    u1, u2 = fzn/fxn, fzn/fwn
    eta = 1/(1 + 2u3*u4^2)/(1 - u2)

    secynzn = secant(fyn, fzn, yn, zn)
    secxnzn = secant(fxn, fzn, xn, zn)
    sec = secynzn - secxnyn + secxnzn
    ((sec == 0.0) | isnan(sec)) && throw(ConvergenceFailed("Division by 0"))

    an = zn - eta * fzn / sec
    fan = f(an)

    u5, u6 = fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 + u2^2*u3 + 3u1*u4^2*(u3^2 - u4^2) / secxnyn
    secynan = secant(fyn, fan, yn, an)
    secznan = secant(fzn, fan, zn, an)
    sec = secynan * secznan
    ((sec == 0.0) | isnan(sec)) && throw(ConvergenceFailed("Division by 0"))
    
    xn1 = zn  - sigma * secynzn * fan / sec

    xn1

end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
## seems faster than order 5 and 8 and more robust than order 2 method
function thukral_update8(f::Function, x0::Real, ftol::Real;
                         beta::Real=1,
                         j::Int=1, k::Int=1, l::Int=1 
                         )

    
    xn = x0
    fxn = f(xn)
    
    abs(fxn) <= ftol && throw(StateConverged(xn))

    ## first step is of size fxn. If that is big, this whole thing has
    ## problems with convergence. Here we replace with a step based on
    ## the central difference
    if abs(fxn/beta) > 1e-1
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end
    

    wn = xn + fxn/beta          # beta can tune how large first step is.
    fwn = f(wn)
    abs(fwn) <= ftol && throw(StateConverged(wn))
    abs(fwn - fxn) == 0.0 && throw(ConvergenceFailed("Division by 0"))

    yn = xn - fxn*fxn/(fwn-fxn)
    fyn = f(yn)
    abs(fyn) <= ftol && throw(StateConverged(yn))
    
    phi = j == 1 ? 1.0 /(1.0 - fyn/fwn) :  (1.0 + fyn/fwn)
    
    abs(fxn - fyn) == 0.0 && throw(ConvergenceFailed("Division by 0"))
    secxnyn = secant(fxn, fyn, xn, yn)
    ((secxnyn == 0.0) | isnan(secxnyn)) && throw(ConvergenceFailed("Division by 0"))
    
    zn = yn - phi * fyn / secxnyn
    fzn = f(zn)

    secynzn = secant(fzn, fyn, zn, yn)
    secxnzn = secant(fzn, fxn, zn, xn)

    rynxn, rynwn, rznwn = fyn/fxn, fyn/fwn, fzn/fwn
    omega = k == 1 ?  1.0 / (1- rznwn) : 1 + rznwn  + rznwn*rznwn
    psi = (l == 1) ? 1 - 2*rynwn*rynwn*rynxn : 1.0 / (1 + 2*rynwn*rynwn*rynxn)

    ## return inc, not new update
    sec = (secynzn - secxnyn + secxnzn)
    ((sec == 0.0) | isnan(sec)) && throw(ConvergenceFailed("Division by 0"))
    
    zn - omega * psi * fzn / sec

end

## Fifth order method
## Numer Algor
## DOI 10.1007/s11075-010-9434-5
## Fifth-order iterative method for finding multiple roots of
##   nonlinear equations
## Xiaowu Li·Chunlai Mu· Jinwen Ma ·Linke Hou
function LiMuMaHou_update(f::Function, xn::Real, ftol::Real; kwargs...)
    
    fxn = f(xn)
    
    ## this is poor if fxn >> 0; we use a hybrid approach with a step
    ## based on f/f' with f' estimated by central difference if fxn is
    ## too big
    if abs(fxn) > 1e-1
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end

    gxn = approx_deriv(f, fxn, xn, ftol)
    gxn == 0.0 && throw(ConvergenceFailed("Division by 0"))
    
    yn = xn - fxn / gxn
    fyn = f(yn)
    abs(fyn) <= ftol && throw(StateConverged(yn))
    
    zn = yn - fyn / gxn
    fzn = f(zn)
    abs(fzn) <= ftol && throw(StateConverged(zn))
    
    sec = F(f, fxn, fyn, fzn, xn, yn, zn, ftol)
    sec == 0.0 && throw(ConvergenceFailed("Division by 0"))

    zn - fzn / sec
end
 
## Some special cases

## Order 2 method
## http://en.wikipedia.org/wiki/Steffensen's_method
function steffensen(g::Function, x0::Real; 
                    ftol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                    xtol::Real =  4.0 * eps(one(eltype(float(x0)))),
                    maxeval::Int=20, 
                    verbose::Bool=false,
                    kwargs...)

    ## we use a fixed point method
    function f(x, ftol) 

        gx = g(x)
        abs(gx) <= ftol && throw(StateConverged(x))

        dg = g(x + gx) - gx
        dg == 0.0 && throw(ConvergenceFailed("Division by 0"))
        
        x - gx * gx / dg
    end

    p0, p = x0, Inf
    verbose && println("p0=$p0, p=$p, ctr=0")
    
    try
        for i in 1:maxeval
            p1 = f(p0, ftol)
            p2 = f(p1, ftol)
            
            ## Aitkens step
            d = (p2 - 2*p1 + p0)
            d == 0.0 && throw(ConvergenceFailed("Division by 0"))

            p=p0 - (p1 - p0)*(p1 - p0) / d

            ## fixed point?
            verbose && println("p0 = $p0, p=$p, ctr=$i")
            abs(p-p0) <= ftol && throw(StateConverged(p))


            p0=p
        end
        throw(ConvergenceFailed("More than $maxeval iterations before convergence"))

    catch e
        if isa(e, StateConverged)
            e.x0
        else
            rethrow(e)
        end
    end
end


"""
Main interface for derivative free methods

* `f` a scalar function `f:R -> R`. Trying to find solution to `f(x) = 0`.

* `x0` initial guess for root. Iterative methods need a reasonable
initial starting point.

* `ftol`. Stop iterating when |f(xn)| <= ftol.

* `xtol`. Stop iterating when |xn+1 - xn| <= xtol.

* `xtolrel`. Stop iterating when |xn+1 - xn| <= |xn| * xtolrel

* `maxeval`. Stop iterating if more than this many steps, throw error.

* `order`. One of 0, 1, 2, 5, 8, or 16. Specifies which algorithm to
use. Default is 0 for the slower, more robust SOLVE function.  Order 8
is faster and pretty robust. The order 2 Steffensen method can be even
faster and use less memory, but is more sensitive to the initial
guess.

* `verbose`. If true, will print out each step taken

* `kwargs...` For order 8, there are some parameters that can be
specified to change the algorithm. In particular, one can specify beta
which controls the size of the first step in an approximate
derivative: (f(x0 + f(x(0))/beta) - f(x0).

The file test/test_fzero2 generates some timed comparisons

Because these are derivative free, they can be used with functions
defined by automatic differentiation
e.g., find critical point of f(x) = x^2
```
fzero(D(x -> x^2), 1)
```
"""
function derivative_free(f::Function, x0::Real;
                         ftol::Real = 10.0 * eps(one(eltype(float(x0)))),
                         xtol::Real =  4.0 * eps(one(eltype(float(x0)))),
                         xtolrel::Real = 10.0 * eps(one(eltype(float(x0)))),
                         maxeval::Int = 30,
                         verbose::Bool=false,
                         order::Int=0, #0, 1, 2, 5, 8 or 16
                         kwargs...      # pass to thukral_update 8, these being beta,j,k,l
                         )

    if order == 16
        update(f::Function, x::Real, ftol) = thukral_update16(f, x, ftol; kwargs...)
    elseif order == 8
        update(f::Function, x::Real, ftol) = thukral_update8(f, x, ftol; kwargs...)
    elseif order == 5
        update(f::Function, x::Real, ftol) = LiMuMaHou_update(f, x, ftol; kwargs...)
    elseif order == 2
        return(steffensen(f, x0, ftol=ftol, xtol=xtol, maxeval=maxeval,
                          verbose=verbose, kwargs...))
    elseif order == 1
        x1 = x0 + min(1e-2, abs(f(x0)))
        return(secant_method(f, x0, x1; ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, verbose=verbose, kwargs...))
    else
        return(SOLVE(f, x0; ftol=ftol, maxeval=maxeval, verbose=verbose, kwargs...))
    end



    xn, xn1 = x0, Inf

    try
        for i in 1:maxeval
            fxn = f(xn)
            abs(fxn) <= ftol && throw(StateConverged(xn))

            xn1 = update(f, xn, ftol)
            isinf(xn1) &&  throw(ConvergenceFailed("sequence diverged"))
            isnan(xn1) &&  throw(ConvergenceFailed("sequence diverged"))
            
            del = abs(xn1 - xn)

            abs(del) <= max(xtol, abs(xn1) * xtolrel) && throw(StateConverged(xn1))

            xn = xn1

            verbose && println("x_$i = $xn;\t f(x_$i) = $(f(xn))")
        end
        throw(ConvergenceFailed("More than $maxeval iterations before convergence"))
    catch e
        if isa(e, StateConverged)
            e.x0
        else
            rethrow(e)
        end
    end
end


