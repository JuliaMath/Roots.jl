## http://www.hindawi.com/journals/ijmms/2012/493456/
## Rajinder Thukral
## very fast (8th order) root finder like newton, but no derivative needed.
##
## Abstract
## A new family of eighth-order derivative-free methods for solving
## nonlinear equations is presented. It is proved that these methods have
## the convergence order of eight. These new methods are derivative-free
## and only use four evaluations of the function per iteration. In fact,
## we have obtained the optimal order of convergence which supports the
## Kung and Traub conjecture. Kung and Traub conjectured that the
## multipoint iteration methods, without memory based on n evaluations
## could achieve optimal convergence order of . Thus, we present new
## derivative-free methods which agree with Kung and Traub conjecture for
## . Numerical comparisons are made to demonstrate the performance of the
## methods presented.
##
## Theorem 2.4. Assume that the function for an open interval D has a
## simple root . Let f: R -> R be sufficiently smooth in the interval D,
## the initial approximation is sufficiently close to then the order of
## convergence of the new derivative-free method defined by (2.7) is
## eight.
##
## 
## some tests
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
## Because this is derivative free, it can be used with functions
## defined by automatic differentiation
## e.g., find minimum of f(x) = x^2
## julia> thukral(D(x -> x^2), 1)

function thukral(f::Function, xo::Real; 
                 tol::Real=1e-14, max_iter::Int=20,
                 verbose::Bool=false,    # return ctr
                 i::Int=1,      # in {1, 2, 3, ...}
                 j::Int=1, k::Int=1, l::Int=1 # in {1,2}
                 )

    function update(f::Function, xn::Real, tol::Real)
        ## i,j,k,l = 1,1 ,1 ,1 ## in {1,2}
        
        beta = 1/i
        fxn = f(xn)

        wn = xn + (1/beta)*fxn
        fwn = f(wn)

        if abs(fwn) < tol return(wn) end

        yn = xn - fxn^2/(fwn - fxn)
        fyn = f(yn)

        if abs(fyn) < tol return(yn) end
        
        phi = ((1 - fyn/fwn)^(-1),  (1 + fyn/fwn))

        zn = yn - phi[j]*( (xn-yn)/(fxn - fyn) ) * fyn
        fzn = f(zn)

        if abs(fzn) < tol return(zn) end
        
        omega = ( (1- fzn/fwn)^(-1), 1 + fzn/fwn  + (fzn/fwn)^2)
        psi = (1 - 2*fyn^3/(fwn^2*fxn), (1 + 2*fyn^3/(fwn^2*fxn))^(-1) )

        zn - omega[k]*psi[l]*((fzn - fyn)/(zn-yn) - (fyn - fxn)/(yn - xn) + (fzn - fxn)/(zn-xn))^(-1)*fzn
    end
    update(f::Function, x::Ad, tol::Real) = update(f, x.val, tol)

    xn = xo
    ctr = 1

    while (abs(f(xn)) > tol) & (ctr < max_iter)
        xn = update(f, xn, tol)
        ctr += 1
    end

    if ctr == max_iter
        throw("Method did not converge in $max_iter")
    else
        verbose ? (xn, ctr) : xn
    end
end
