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
    abs(a-b) == 0.0 && throw(StateConverged(b))
    (fa - fb)/(a-b)
end

function approx_deriv(f::Function, fx::Real, x::Real)
    fx == 0.0 && throw(StateConverged(x))
    (f(x + fx) - fx)/fx
end

function secant2(f, fz, fx, z, x) 
    abs(x-z) == 0.0 && throw(StateConverged(x))
    (approx_deriv(f, fx, x) - secant(fz, fx, z, x))/(x-z)
end

function F(f::Function, fx, fy, fz, x, y, z)
    secant(fz, fy, z, y)  + secant2(f, fz, fx, z, x) * (z-y)
end

isissue(x) = isnan(x) || isinf(x) || x == 0.0 # check approx derivatives with this
isissue(x, y) = isissue(x) || abs(x/y) > 1e1 || abs(x/y) < 1e-1

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
function thukral_update16(f::Function, x0::Real; kwargs...)

    xn = x0
    fxn = f(xn)

    wn = xn + fxn
    fwn = f(wn)

    fxn == 0.0 && throw(StateConverged(xn))
    fwn == 0.0 && throw(StateConverged(wn))

    ## first step is of size fxn. If that is big, this whole thing has
    ## problems with convergence. Here we replace with a step based on
    ## the central difference
    if abs(fxn) > 1e-1 
        h = 1e-6
        return(secant_step(f(xn+h), f(xn-h), xn+h, xn-h))
    end

    secxnwn = secant(fxn, fwn, xn, wn)
    
    yn = xn - fxn/secxnwn
    fyn = f(yn)

    secxnyn = secant(fxn, fyn, xn, yn)
    secwnyn = secant(fwn, fyn, wn, yn)

    u3, u4 = fyn/fxn, fyn/fwn
    phi1, phi2, phi3 = 1/(1 - u4), 1 + u4 + u4^2, secxnwn / secwnyn


    zn = yn - phi3 * fyn / secxnyn
    fzn = f(zn)

    u1, u2 = fzn/fxn, fzn/fwn
    eta = 1/(1 + 2u3*u4^2)/(1 - u2)

    secynzn = secant(fyn, fzn, yn, zn)
    secxnzn = secant(fxn, fzn, xn, zn)
    sec = secynzn - secxnyn + secxnzn
    sec == 0.0 && throw(StateConverged(zn))

    an = zn - eta * fzn / sec
    fan = f(an)

    u5, u6 = fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 + u2^2*u3 + 3u1*u4^2*(u3^2 - u4^2) / secxnyn
    secynan = secant(fyn, fan, yn, an)
    secznan = secant(fzn, fan, zn, an)

    xn1 = zn  - sigma * secynzn * fan / secynan / secznan

    xn1

end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
## seems faster than order 5 and 8 and more robust than order 2 method
function thukral_update8(f::Function, x0::Real;
                         beta::Real=1,
                         j::Int=1, k::Int=1, l::Int=1 
                         )

    
    xn = x0
    fxn = f(xn)
    
    fxn == 0.0 && throw(StateConverged(xn))

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
    fwn == 0.0 && throw(StateConverged(wn))
    abs(fwn - fxn) == 0.0 && throw(StateConverged(wn))

    yn = xn - fxn*fxn/(fwn-fxn)
    fyn = f(yn)

    phi = j == 1 ? 1.0 /(1.0 - fyn/fwn) :  (1.0 + fyn/fwn)
    
    abs(fxn - fyn) == 0.0 && throw(StateConverged(yn))
    secxnyn = secant(fxn, fyn, xn, yn)
    zn = yn - phi * fyn / secxnyn
    fzn = f(zn)


    secynzn = secant(fzn, fyn, zn, yn)
    secxnzn = secant(fzn, fxn, zn, xn)


    rynxn, rynwn, rznwn = fyn/fxn, fyn/fwn, fzn/fwn
    omega = k == 1 ?  1.0 / (1- rznwn) : 1 + rznwn  + rznwn*rznwn
    psi = (l == 1) ? 1 - 2*rynwn*rynwn*rynxn : 1.0 / (1 + 2*rynwn*rynwn*rynxn)

    ## return inc, not new update
    sec = (secynzn - secxnyn + secxnzn)
    sec == 0.0 && throw(StateConverged(zn))
    zn - omega * psi * fzn / sec

end

## Fifth order method
## Numer Algor
## DOI 10.1007/s11075-010-9434-5
## Fifth-order iterative method for finding multiple roots of
##   nonlinear equations
## Xiaowu Li·Chunlai Mu· Jinwen Ma ·Linke Hou
function LiMuMaHou_update(f::Function, xn::Real; kwargs...)
    
    fxn = f(xn)

    ## this is poor if fxn >> 0; we use a hybrid approach with a step
    ## based on f/f' with f' estimated by central difference if fxn is
    ## too big
    if abs(fxn) > 1e-1
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end

    gxn = approx_deriv(f, fxn, xn)
    yn = xn - fxn / gxn
    fyn = f(yn)

    zn = yn - fyn / gxn
    fzn = f(zn)

    zn - fzn / F(f, fxn, fyn, fzn, xn, yn, zn)
end
 
## Some special cases

## Order 2 method
## http://en.wikipedia.org/wiki/Steffensen's_method
function steffensen(g::Function, x0::Real; 
                    tol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                    delta::Real =  4.0 * eps(one(eltype(float(x0)))),
                    max_iter::Int=100, 
                    verbose::Bool=false,
                    kwargs...)

    ## we use a fixed point method
    function f(x) 

        gx = g(x)
        gx == 0.0 && throw(StateConverged(x))

        gp = (g(x + gx) - gx) / gx
        gp == 0.0 && throw(StateConverged(x))

        x - gx/gp
    end

    p0, p = x0, Inf

    try
        for i in 1:max_iter
            verbose && println("p0 = $p0, p=$p, ctr=$(i-1)")

            p1 = f(p0)

            abs(p1) <= tol && throw(StateConverged(p0))
            isnan(p1) && throw(StateConverged(p0))

            p2 = f(p1)
            isnan(p2) && throw(StateConverged(p1))

            ## Aitkens step
            p=p0-(p1-p0)^2/(p2-2*p1+p0)

            abs(p-p0) <= tol && throw(StateConverged(p))
            p0=p
        end
        throw(ConvergenceFailed())

    catch e
        if isa(e, StateConverged)
            e.x0
        else
            throw(ConvergenceFailed())
        end
    end
end

function steffensen_update(f::Function, x::Real; kwargs...)
    fx = f(x)
    if abs(fx) > 1e-1 
        h = 1e-6
        return(secant_step(f(x+h), f(x-h), x+h, x-h))
    end

    den = f(x+fx) - fx
    den == 0.0 && throw(StateConverged(x))

    x - fx*fx/den
end

## Order 1 secant method
function secant_method(f::Function, x0::Real, x1::Real;
                tol::Real   = 10.0 * eps(one(eltype(float(x1)))),
                delta::Real =  zero(x1),
                max_iter::Int=100, 
                verbose::Bool=false,
                kwargs...)

    a, b, fa, fb = x0, x1, f(x0), f(x1)
    
    try
        fb == 0 && throw(StateConverged(b))

        for i in 1:max_iter
            verbose && println("a=$a, b=$b, ctr=$(i-1)")

            inc = fb/secant(fa, fb, a, b)
            a, b = b, b-inc
            fa, fb = fb, f(b)
            abs(inc) <= delta && throw(StateConverged(b))
            abs(fb) <= tol && throw(StateConverged(b))

        end

        throw(ConvergenceFailed())

    catch e
        if isa(e, StateConverged)
            e.x0
        else
            throw(e)
        end
    end
end

secant_step(fa, fb, a, b) = b - fb / secant(fa, fb, a, b)

## Main interface
##
## @param f a scalar function f:R -> R. Trying to find solution to f(x) = 0.
## @param x0 initial guess for root. Iterative methods need a
##        reasonable initial starting point.
## @param tol. Stop iterating when |f(xn)| <= tol.
## @param delta. Stop iterating when |xn+1 - xn| <= delta.
## @param max_iter. Stop iterating if more than this many steps, throw error.
## @param order. One of 2, 5, 8, or 16. Specifies which algorithm to
##        use. Default is 2 (Steffensen's method), as this is generally faster, but
##        it is less robust to the initial guess
## @param verbose. If true, will print out each step taken
## @param kwargs... For order 8, there are some parameters that can be
##        specified to change the algorithm. In particular, one can specify
##        beta which controls the size of the first step in an approximate
##        derivative: (f(x0 + f(x(0))/beta) - f(x0).
## 
##
## The file test/test_fzero2 generates some timed comparisons
##
## Because these are derivative free, they can be used with functions
## defined by automatic differentiation
## e.g., find critical point of f(x) = x^2
## julia> fzero(D(x -> x^2), 1)


function derivative_free(f::Function, x0::Real;
                 tol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                 delta::Real =  4.0 * eps(one(eltype(float(x0)))),
                 max_iter::Int = 200,
                 verbose::Bool=false,
                 order::Int=2, # 2, 5, 8 or 16
                 kwargs...      # pass to thukral_update 8, these being beta,j,k,l
                 )

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x; kwargs...)
    elseif order == 5
        update(f::Function, x::Real) = LiMuMaHou_update(f, x; kwargs...)
    elseif order == 2
#        update(f::Function, x::Real) = steffensen_update(f, x)
        ## much faster
        return(steffensen(f, x0, tol=tol, delta=delta, max_iter=max_iter,
                          verbose=verbose, kwargs...))
    else
        update(f::Function, x::Real) = thukral_update8(f, x; kwargs...)
    end



    xn, xn1, del = x0, Inf, Inf

    try
        for i in 1:max_iter
            abs(f(xn)) <= tol && throw(StateConverged(xn))
            abs(del) <= delta && throw(StateConverged(xn))

            xn1 = update(f, xn)
            del = abs(xn1 - xn)
            xn = xn1
        
            verbose && println("x_$i = $xn;\t f(x_$i) = $(f(xn))")
        end
        throw(ConvergenceFailed())
    catch e
        if isa(e, StateConverged)
            e.x0
        else
            throw(e)
        end
    end
end

## bracket answer, if algorithm leaves bracket, then switch to
## bracketing algorithm
function derivative_free_bracket(f::Function, x0::Real, bracket::Vector;
                         tol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                         delta::Real =  zero(x0),
                         max_iter::Int = 200,
                         verbose::Bool=false,
                         order::Int=8,
                         kwargs...
                         )
    

    a, b = sort(bracket[1:2])

    f(a) * f(b) < 0 || error("The bracket does not bracket a root")
    a <= x0 <= b || error("x0 not in [a,b]")

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x; kwargs...)
    elseif order == 5
        update(f::Function, x::Real) = LiMuMaHou_update(f, x; kwargs...)
    elseif order == 2
        update(f::Function, x::Real) = steffensen_update(f, x; kwargs...)
    else
        update(f::Function, x::Real) = thukral_update8(f, x; kwargs...)
    end

    try
        xn = x0
        del = Inf
        for i in 1:max_iter
            xn1 = update(f, xn)
        
            if isissue(xn1)   # no valid step returned
                xn1 = a + (b-a)/2
            end
            if xn1 < a || xn1 > b   # outside of bracket 
                xn1 = a + (b-a)/2
            end
            del = xn1 - xn
            xn = xn1
            ## check tolerances
            abs(f(xn)) <= tol && throw(StateConverged(xn))
            abs(del) <= delta && throw(StateConverged(xn))
            
            ## update bracket
            f(a) * f(xn) < 0 ? (b = xn) : (a=xn)
            
            verbose && println("xn1=", xn, " step=", i) 
        end
        throw(ConvergenceFailed())

    catch e
        if isa(e, StateConverged)
            e.x0
        else
            throw(e)
        end
    end
end  

