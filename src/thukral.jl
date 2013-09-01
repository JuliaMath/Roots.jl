## Various higher order, derivative free root finding algorithms
## These are all basically based on an update step which is 3 or 4 newton's method calls
## of the type xn1 = xn - w f(xn)/f'(xn) with different choices for the weights
## and different derivative-free approximations for f'(xn).
##
## The 8th order thukral method is consistently faster and uses fewer allocated bytes

## some helpers
secant(fa::Real, fb::Real, a::Real, b::Real) = (fa - fb)/(a-b)
secant(fa::Ad, fb::Ad, a::Ad, b::Ad) = (fa.val - fb.val)/(a.val-b.val)
function secant(fa::Ad, fb, a, b) 
    fa, fb, a, b = promote(fa, fb, a, b)
    secant(fa, fb, a, b)
end
approx_deriv(f::Function, fx::Real, x::Real) = (f(x + fx) - fx)/fx
function secant2(f, fz, fx, z, x) 
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
function thukral_update16(f::Function, x0::Real, delta::Real; kwargs...)

    xn = x0
    fxn = f(xn)

    wn = xn + fxn
    fwn = f(wn)

    ## first step is of size fxn. If that is big, this whole thing has problems
    ## with convergence. Here we replace with a step based on the approximate derivative
    if abs(fxn) > 1e-2
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end

    adiff1 = secxnwn = secant(fxn, fwn, xn, wn)
    isissue(adiff1) && return(xn) # can't improve
    
    inc = fxn / adiff1
    isnan(inc) && return(xn)

    yn = xn - fxn/adiff1
    abs(inc) < delta && return(yn)

    fyn = f(yn)

    adiff2 = secxnyn = secant(fxn, fyn, xn, yn)
    secxnyn = secant(fxn, fyn, xn, yn)
    secwnyn = secant(fwn, fyn, wn, yn)

    u3, u4 = fyn/fxn, fyn/fwn
    phi1, phi2, phi3 = 1/(1 - u4), 1 + u4 + u4^2, secxnwn / secwnyn

    inc = phi3 * fyn / (!isissue(adiff2) ? adiff2 : adiff1)
    isnan(inc) && return(yn)

    zn = yn - inc
    abs(inc) < delta && return(zn)
    fzn = f(zn)

    u1, u2 = fzn/fxn, fzn/fwn

    eta = 1/(1 + 2u3*u4^2)/(1 - u2)

    secynzn = secant(fyn, fzn, yn, zn)
    secxnzn = secant(fxn, fzn, xn, zn)
    adiff3 = secynzn - secxnyn + secxnzn

    an = zn - eta * fzn / (!isissue(adiff3) ? adiff3 : (!isissue(adiff2) ? adiff2 : adiff1))
    fan = f(an)

    u5, u6 = fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 + u2^2*u3 + 3u1*u4^2*(u3^2 - u4^2) / secxnyn
    secynan = secant(fyn, fan, yn, an)
    secznan = secant(fzn, fan, zn, an)
    adiff4 = secynzn / secynan / secznan

    xn1 = zn - sigma * fan / (!isissue(adiff4) ? adiff4 :(!isissue(adiff3) ? adiff3 : (!isissue(adiff2) ? adiff2 : adiff1)))

    xn1

end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
## We use this as the default. Seems faster than update16 and takes less memory
function thukral_update8(f::Function, x0::Real, delta::Real;
                         beta::Real=1,
                         j::Int=1, k::Int=1, l::Int=1 
                         )

    
    xn = x0
    fxn = f(xn)
    

    ## this is poor if fxn >> 0; we use a hybrid approach with
    ## a step based on f/f' with f' estimated by central difference if fxn is too big
    if abs(fxn/beta) > 1e-2
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end
    

    wn = xn + fxn/beta          # beta can tune how large first step is.
    fwn = f(wn)
    adiff1 = (fwn - fxn) / (fxn/beta) 
    inc = (1 / beta) * fxn/adiff1
    
    isissue(adiff1) && return(xn)

    yn = xn - inc
    abs(inc) < delta && return(yn)
    
    fyn = f(yn)
    secxnyn = secant(fxn, fyn, xn, yn)
    adiff2 = secxnyn
    
    phi = j == 1 ? 1.0 /(1.0 - fyn/fwn) :  (1.0 + fyn/fwn)
    
    inc = phi * fyn / (isissue(adiff2, adiff1)  ? adiff1 : adiff2)
    isnan(inc) && return(yn)

    zn = yn - inc
    abs(inc) < delta && return(zn)
    fzn = f(zn)

    secynzn = secant(fzn, fyn, zn, yn)
    secxnzn = secant(fzn, fxn, zn, xn)
    rynxn, rynwn, rznwn = fyn/fxn, fyn/fwn, fzn/fwn

    omega = k == 1 ?  1.0 / (1- rznwn) : 1 + rznwn  + rznwn*rznwn
    psi = (l == 1) ? 1 - 2*rynwn*rynwn*rynxn : 1.0 / (1 + 2*rynwn*rynwn*rynxn)

    adiff3 = secynzn - secxnyn + secxnzn

    inc = omega * psi * fzn / (isissue(adiff3, adiff2) ? (isissue(adiff2) ? adiff1 : adiff2) : adiff3)

    ##
    zn - inc

end

## Fifth order method
## Numer Algor
## DOI 10.1007/s11075-010-9434-5
## Fifth-order iterative method for finding multiple roots of nonlinear equations
## Xiaowu Li·Chunlai Mu· Jinwen Ma ·Linke Hou
function LiMuMaHou_update(f::Function, xn::Real, delta::Real; kwargs...)
    
    fxn = f(xn)
    gxn = approx_deriv(f, fxn, xn)

    ## this is poor if fxn >> 0; we use a hybrid approach with
    ## a step based on f/f' with f' estimated by central difference if fxn is too big
    if abs(fxn) > 1e-1
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end
    

    inc = fxn / gxn
    isnan(inc) && return(xn)

    yn = xn - inc
    abs(inc) < delta && return(yn)
    

    fyn = f(yn)
    inc = fyn / gxn
    isnan(inc) && return(yn)

    zn = yn - inc



    fzn = f(zn)
    inc = fzn / F(f, fxn, fyn, fzn, xn, yn, zn)
    
    isnan(inc) && return(zn)

    zn - inc
end


## Main interface
##
## @param f a scalar function f:R -> R. Trying to find solution to f(x) = 0.
## @param x0 initial guess for root. Iterative methods need a reasonable initial starting point.
## @param tol. Stop iterating when |f(xn)| <= tol.
## @param delta. Stop iterating when |xn+1 - xn| <= delta.
## @param max_iter. Stop iterating if more than this many steps, throw error.
## @param order. One of 5, 8, or 16. Specifies which algorithm to use. Default is 8 which seems to win both in 
##        speed of execution and in memory consumption.
## @param verbose. If true, will print out each step taken
## @param kwargs... For order 8, there are some parameters that can be
##        specified to change the algorithm. In particular, one can specify
##        beta which controls the size of the first step in an approximate
##        derivative: (f(x0 + f(x(0))/beta) - f(x0).
## 
## We have 5, 8 and 16 order methods. Empirically it seems 16 converges sometimes when 8 does not, though 8 is a bit faster.
## 
## some tests. (See also http://ir.igsnrr.ac.cn/bitstream/311030/8840/1/%E4%BE%AF%E9%BA%9F%E7%A7%91(SCI)2.pdf)
## ------
## julia> thukral(x -> (x-2)*(x^10 + x + 1)*exp(-x-1), 1.9)
## (2.0,4)
##
## julia> thukral(x -> x^11 + x + 1, -1)
## (-0.844397528792023,4)
##
## julia> thukral(u -> sin(u)^2 - u^2 + 1, 1)
## (1.4044916482153413,3)
##
## ## some comparison for a tricky function
## julia> f(x) = (log(x) + sqrt(x^4+1)-2)^7
## julia> newton(f, 1, verbose=true)            ## 26 steps; 0.064031741 seconds (148824 bytes allocated)
## julia> thukral(f, 1, order=16, verbose=true) ## 15 steps; 0.000126401 seconds (23320 bytes allocated)
## julia> thukral(f, 1, order=8, verbose=true)  ## 18 steps; 8.7874e-5 seconds (13704 bytes allocated)
## julia> thukral(f, 1, order=5, verbose=true)  ## 11 steps; 9.7764e-5 seconds (9992 bytes allocated)
## julia> fzero(f, [1,2])                       ## 42 steps; 0.000698867 seconds (159024 bytes allocated)
##
## Can often get more accuracy with relatively little cost by using BigFloat:
## julia> f(x) = (8x*exp(-x^2) -2x - 3)^8
## julia> @time fzero(f, x0) - -1.7903531791589544 ## only 3e-3! High multiplicity, |f(xstar)| < 1e-15
## julia> @time fzero(f, BigFloat(x0)) - -1.7903531791589544 ## now 6e-11. Took 100 times as long
## julia> @time newton(f, x0) |> f ## about same time as with BigFloat and similar accuracy; 27 steps
##
## Because these are derivative free, they can be used with functions
## defined by automatic differentiation
## e.g., find minimum of f(x) = x^2
## julia> thukral(D(x -> x^2), 1)
##
## ## This won't always be perfect though. For example, the function
## ## A(t) = 100sin(t)*cos(t) + 100*sin(t) will only get to within
## ## 1e-14 when trying to find fzero(D(A), 1):
## julia> A(t) = 100sin(t)*cos(t) + 100*sin(t) 
## julia> fzero(D(A), 1) |> D(A)
## 2.842170943040401e-14

function thukral(f::Function, x0::Real;
                 tol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                 delta::Real =  4.0 * eps(one(eltype(float(x0)))),
                 max_iter::Int = 100,
                 verbose::Bool=false,
                 order::Int=8, # 5, 8 or 16
                 kwargs...      # pass to thukral_update 8, these being beta,j,k,l
                 )

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x, delta; kwargs...)
    elseif order == 5
        update(f::Function, x::Real) = LiMuMaHou_update(f, x, delta; kwargs...)
    else
        update(f::Function, x::Real) = thukral_update8(f, x, delta; kwargs...)
    end
    update(f::Function, x::Ad) = update(f, x.val)


    xn, xn1, del = x0, Inf, Inf
    cvg = false

    for i in 1:max_iter
        if abs(f(xn)) <= tol
            cvg = true
            break
        end
        if abs(del) <= delta
            cvg = true
            break
        end

        xn1 = update(f, xn)
        del = abs(xn1 - xn)
        xn = xn1
        
        verbose && println("xn= $xn, f(xn) = $(f(xn)), del=$del, ctr=$i")
    end

    cvg || error("Did not converge in $max_iter steps")

    xn
end

## bracket answer, if algorithm leaves bracket, then switch to bracketing algorithm
function thukral_bracket(f::Function, x0::Real, bracket::Vector;
                         tol::Real   = 10.0 * eps(one(eltype(float(x0)))),
                         delta::Real =  4.0 * eps(one(eltype(float(x0)))),
                         max_iter::Int = 100,
                         verbose::Bool=false,
                         order::Int=8,
                         kwargs...
                         )
    

    a, b = sort(bracket[1:2])
    f(a) * f(b) < 0 || error("The bracket does not bracket a root")
    a <= x0 <= b || error("x0 not in [a,b]")

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x, delta; kwargs...)
    elseif order == 5
        update(f::Function, x::Real) = LiMuMaHou_update(f, x, delta; kwargs...)
    else
        update(f::Function, x::Real) = thukral_update8(f, x, delta; kwargs...)
    end
    update(f::Function, x::Ad) = update(f, x.val)

    cvg = false
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
        if abs(f(xn)) <= tol || abs(del) < delta
            cvg=true
            break
        end
        
        ## update bracket
        f(a) * f(xn) < 0 ? (b = xn) : (a=xn)

        verbose ? println("xn1=", xn, " step=", i) : nothing
    end

    cvg || error("Method did not converge in $max_iter")
    
    xn
    
end  
