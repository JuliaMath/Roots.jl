## Various higher order, derivative free root finding algorithms

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
function thukral_update16(f::Function, x0::Real, tol::Real; beta::Real=1)
    xn = x0
    fxn = f(xn)

    ## can have issues if fxn is too big, we replace here with hybrid step
    if abs(fxn) > 1e-2 * beta
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end

    wn = xn + fxn / beta
    fwn = f(wn)

    inc = fxn / secant(fwn, fxn, wn, xn)
    abs(inc) <= tol && return(xn)
    

    yn = xn - inc
    fyn = f(yn)


    phi3 = secant(fxn, fwn, xn, wn) / secant(fyn, fwn, yn, wn)
    inc = phi3 * fyn / secant(fxn, fyn, xn , yn)
    abs(inc) <= tol && return(yn)

    zn = yn - inc
    fzn = f(zn)


    u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn 
    eta = 1.0 / (1 + 2*u3*u4*u4) / (1 - u2)

    inc = eta * fzn / (secant(fyn, fzn, yn, zn) - secant(fxn, fyn, xn, yn) + secant(fxn, fzn, xn, zn))
    abs(inc) <= tol && return(zn)

    an = zn -  inc
    fan = f(an)

    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1 *u3*u4*u4 + u5 + u6 + u1*u1*u4 + u2*u2*u3 + 3*u1*u4*u4*(u3*u3 - u4*u4) / secant(fxn, fyn, xn, yn)
    
    inc = sigma * secant(fyn, fzn, yn, zn) * fan / secant(fyn, fan, yn, an) / secant(fzn, fan, zn, an)
    zn - inc
end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) derivative free iterative root finder.
## We use this as the default. Seems faster than update16 and takes less memory
function thukral_update8(f::Function, x0::Real, tol::Real;
                         beta::Real=1,
                         j::Int=1, k::Int=1, l::Int=1 
                         )
        
    xn = x0
    fxn = f(xn)

    ## first step is based on f'(xn) approx f(xn + f(xn))/f(xn)
    ## can have issues if fxn is too big, we replace here with hybrid step
    if  abs(fxn) > 1e-2 * beta
        h = 1e-6
        fp = (f(xn+h) - f(xn-h)) / (2h)
        return (xn - fxn/fp)
    end

    wn = xn + fxn/beta
    fwn = f(wn)

    inc =  fxn * fxn / (fwn - fxn)
    abs(inc) <= tol && return(wn)

    yn = xn - inc
    fyn = f(yn)
    
    phi = j == 1 ? 1.0 /(1.0 - fyn/fwn) :  (1.0 + fyn/fwn)

    inc = phi* fyn / secant(fxn, fyn, xn, yn)
    abs(inc) <= tol && return(yn)

    zn = yn - inc
    fzn = f(zn)

    omega = k == 1 ?  1.0 / (1- fzn/fwn) : 1 + fzn/fwn  + fzn*fzn/fwn/fwn
    psi = (l == 1) ? 1 - 2*fyn*fyn*fyn/(fwn*fwn*fxn) : 1.0 / (1 + 2*fyn*fyn*fyn/(fwn*fwn*fxn))

    inc = omega * psi * fzn / (secant(fzn, fyn, zn, yn) - secant(fyn, fxn, yn, xn) + secant(fzn, fxn, zn, xn)) 

    zn - inc
end

## Fifth order method
## Numer Algor
## DOI 10.1007/s11075-010-9434-5
## Fifth-order iterative method for finding multiple roots of nonlinear equations
## Xiaowu Li·Chunlai Mu· Jinwen Ma ·Linke Hou
function LiMuHou_update(f::Function, xn::Real, delta::Real; beta::Real=1.0)
    
    fxn = f(xn)
    gxn = approx_deriv(f, fxn, xn)

    ## this is poor if fxn >> 0; we use a hybrid approach with
    ## a step based on f/f' with f' estimated by central difference if fxn is to big
    if abs(fxn) < 1e-2 
        ## regularly scheduled program
        inc = fxn / gxn
    else
        ## use secant line with h=1e-6?
        h = 1e-6
        fpxn = (f(xn+h) - f(xn-h)) / 2h
        inc = fxn/fpxn
    end
    yn = xn - inc

    fyn = f(yn)
    inc = fyn / gxn
    zn = yn - inc
    abs(inc) <= delta && return(zn)

    fzn = f(zn)
    inc = fzn / F(f, fxn, fyn, fzn, xn, yn, zn)

    zn - inc
end


## Main interface
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
## julia> newton(f, 1, verbose=true) ## 26 steps! elapsed time: 0.064031741 seconds (148824 bytes allocated)
## julia> thukral(f, 1, order-16, verbose=true) ## 12 steps, elapsed time: 0.000199814 seconds (30296 bytes allocated)
## julia> thukral(f, 1, order=8, verbose=true) ## 7 steps, elapsed time: 0.000157921 seconds (14592 bytes allocated)
## julia> fzero(f, [1,2]) ## 42 steps, elapsed time: 0.001219398 seconds (166808 bytes allocated)
##
## Because these are derivative free, they can be used with functions
## defined by automatic differentiation
## e.g., find minimum of f(x) = x^2
## julia> thukral(D(x -> x^2), 1)
##
function thukral(f::Function, x0::Real;
                 tol::Real = 10.0 * eps(1.0),
                 delta::Real = 4 * eps(1.0),
                 max_iter::Int = 100,
                 verbose::Bool=false,
                 order::Int=16, # 5, 8 or 16
                 kwargs...      # pass to thukral_update 8, these being beta,j,k,l
                 )

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x, delta; kwargs...)
    elseif order == 5
        update(f::Function, x::Real) = LiMuHou_update(f, x, delta; kwargs...)
    else
        update(f::Function, x::Real) = thukral_update8(f, x, delta; kwargs...)
    end
    update(f::Function, x::Ad) = update(f, x.val)

    xn, xn1, del = x0, Inf, Inf
    cvg = false

    for i in 0:max_iter
        if abs(f(xn)) <= tol || del <= delta
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
                         delta::Real   = 10.0 * eps(1.0),
                         tol::Real     = 10.0 * eps(1.0),
                         max_iter::Int = 100,
                         verbose::Bool=false,
                         order::Int=16,
                         kwargs...
                         )
    

    a, b = sort(bracket[1:2])
    f(a) * f(b) < 0 || error("The bracket does not bracket a root")
    a <= x0 <= b || error("x0 not in [a,b]")

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x, delta)
    elseif order == 5
        update(f::Function, x::Real) = LiMuHou_update(f, x, delta)
    else
        update(f::Function, x::Real) = thukral_update8(f, x, delta; kwargs...)
    end
    update(f::Function, x::Ad) = update(f, x.val)

    cvg = false
    xn = x0
    for i in 1:max_iter
        xn1 = update(f, xn)
        if isnan(xn1) || isinf(xn1) || xn1 < a || xn1 > b
            xn1 = newton_quadratic(f::Function, a, b, xn, 3)
        end
        xn = xn1

        if abs(f(xn)) <= tol || b - a < delta
            cvg=true
            break
        end
        
        f(a) * f(xn) < 0 ? (b = b - (b-xn)/2) : (a=a + (xn-a)/2)

        verbose ? println("xn=", xn, " step=", i) : nothing
    end

    cvg || throw("Method did not converge in $max_iter")
    
    xn
end  