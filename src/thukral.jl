## 16th order, derivative free root finding algorithm
## http://article.sapub.org/10.5923.j.ajcam.20120203.08.html
## American Journal of Computational and Applied Mathematics
## p-ISSN: 2165-8935    e-ISSN: 2165-8943
## 2012;  2(3): 112-118
## doi: 10.5923/j.ajcam.20120203.08
## New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations
## R. Thukral
## Research Centre, 39 Deanswood Hill, Leeds, West Yorkshire, LS17 5JS, England

secant(fa, fb, a, b) = (fa - fb)/(a-b)
function thukral_update16(f::Function, x0::Real)
    xn = x0
    fxn = f(xn)

    wn = xn + fxn
    fwn = f(wn)
    wn == xn && return(wn)

    yn = xn - fxn / secant(fwn, fxn, wn, xn)
    fyn = f(yn)
    yn == wn && return(yn)

    phi3 = secant(fxn, fwn, xn, wn) / secant(fyn, fwn, yn, wn)

    zn = yn - phi3 * fyn / secant(fxn, fyn, xn , yn)
    fzn = f(zn)
    zn == yn && return(zn)

    u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn 
    eta = 1.0 / (1 + 2*u3*u4*u4) / (1 - u2)

    an = zn - eta * fzn / (secant(fyn, fzn, yn, zn) - secant(fxn, fyn, xn, yn) + secant(fxn, fzn, xn, zn))
    fan = f(an)
    an == zn && return(an)

    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
    sigma = 1 + u1*u2 - u1 *u3*u4*u4 + u5 + u6 + u1*u1*u4 + u2*u2*u3 + 3*u1*u4*u4*(u3*u3 - u4*u4) / secant(fxn, fyn, xn, yn)
    
    zn - sigma * secant(fyn, fzn, yn, zn) * fan / secant(fyn, fan, yn, an) / secant(fzn, fan, zn, an)
end


## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) root finder like newton, but no derivative needed.
function thukral_update8(f::Function, xn::Real, tol::Real;
                         i::Int=1, j::Int=1,
                         k::Int=1, l::Int=1)
        
    beta = 1/i
    fxn = f(xn)
    
    wn = xn + (1/beta)*fxn
    fwn = f(wn)
    
    if abs(fwn) < tol return(wn) end
    
    yn = xn - fxn^2/(fwn - fxn)
    fyn = f(yn)
    
    if abs(fyn) < tol return(yn) end
    
    phi = j == 1 ? (1 - fyn/fwn)^(-1) :  (1 + fyn/fwn)
    
    zn = yn - phi*( (xn-yn)/(fxn - fyn) ) * fyn
    fzn = f(zn)
    
    if abs(fzn) < tol return(zn) end
    
    omega = k == 1 ?  (1- fzn/fwn)^(-1) : 1 + fzn/fwn  + (fzn/fwn)^2
    psi = (l == 1) ? 1 - 2*fyn^3/(fwn^2*fxn) : (1 + 2*fyn^3/(fwn^2*fxn))^(-1) 
    
    zn - omega*psi*((fzn - fyn)/(zn-yn) - (fyn - fxn)/(yn - xn) + (fzn - fxn)/(zn-xn))^(-1)*fzn
end

## Main interface
## 
## We have both 8 and 16, though not sure that 16 performs better than 8.
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
## Because this is derivative free, it can be used with functions
## defined by automatic differentiation
## e.g., find minimum of f(x) = x^2
## julia> thukral(D(x -> x^2), 1)
##
function thukral(f::Function, x0::Real;
                 tol::Real = 4*eps(1.0),
                 delta::Real = 4*eps(1.),
                 max_iter::Int = 100,
                 verbose::Bool=false,
                 order::Int=8, # 8 or 16
                 kwargs...      # pass to thukral_update 8, these being i,j,k,l
                 )

    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x)
        update(f::Function, x::Ad) = thukral_update16(f, x.val)
    else
        update(f::Function, x::Real) = thukral_update8(f, x, tol; kwargs...)
        update(f::Function, x::Ad) = thukral_update8(f, x.val, tol; kwargs...)
    end

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
        
        verbose && println("xn= $xn, f(xn) = $(f(xn)), ctr=$i")
    end

    cvg || error("Did not converge in $max_iter steps")

    xn
end

## bracket answer, if algorithm leaves bracket, then switch to bracketing algorithm
function thukral_bracket(f::Function, x0::Real, bracket::Vector;
                         delta::Real   = 100.0 * eps(1.0),
                         tol::Real     = 100.0 * eps(1.0),
                         max_iter::Int = 100,
                         verbose::Bool=false,
                         order::Int=8,
                         kwargs...
                         )
    

    a, b = sort(bracket[1:2])
    f(a) * f(b) < 0 || error("The bracket does not bracket a root")
    a <= x0 <= b || error("x0 not in [a,b]")


    if order == 16
        update(f::Function, x::Real) = thukral_update16(f, x)
        update(f::Function, x::Ad) = thukral_update16(f, x.val)
    else
        update(f::Function, x::Real, tol::Real) = thukral_update8(f, x, tol; kwargs...)
        update(f::Function, x::Ad, tol::Real) = thukral_update8(f, x.val, tol; kwargs...)
    end

    cvg = false
    xn = x0
    for i in 1:max_iter
        xn1 = update(f, xn, tol)
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
    
end  