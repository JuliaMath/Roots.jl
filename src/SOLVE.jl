
type PossibleExtremaReached
    x0::Real
end

type Dithering
    x0::Real
end

"""

SOLVE button from the HP-34C
follows roughly algorithm described http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf
though some modifications were made.

Goal is to return a value `x` with either:
* `f(x) == 0.0` or 
* `f(prevfloat(x)) * f(nextfloat(x)) < 0`.
if a bracket is found that can be done. 

Otherwise, values with `f(x) <= ftol` are given status of answers
and others are marked as a possible extrema.

If the algorithm gets stuck, a Dithering error is thrown.

While this is slower than fzero with order 2, 5, 8, or 16 it is 
more robust to the initial condition.

"""
function SOLVE(f, x0::Real;
               ftol::Real=10*eps(one(float(x0))),
               xtol::Real=zero(x0),
               xtolrel::Real=zero(x0),
               maxeval::Int=30,
               verbose::Bool=false)

    x0 =  float(x0)
    fx0 = f(x0)
    abs(fx0) <= ftol && return(x0)
    
    try
        ## get initial guess. Use Steffensen guess if reasonable size
        stepsize = max(1/100, min(abs(fx0), abs(x0/100)))
        if fx0 < 0
            g(x) = -f(x)
            secant_method_no_bracket(g, x0 - stepsize, x0, ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, verbose=verbose)
        else
            secant_method_no_bracket(f, x0 - stepsize, x0, ftol=ftol, xtol=xtol, xtolrel=xtolrel, maxeval=maxeval, verbose=verbose)
        end
    catch ex
        if isa(ex, StateConverged)
            return ex.x0
        elseif isa(ex, PossibleExtremaReached)
            guess = ex.x0
            if f(nextfloat(guess)) * f(prevfloat(guess)) <= 0
                return guess
            end
            # heuristic for case where it thinks it is at a minimum
            if abs(f(guess)) < ftol
                return guess
            else
                rethrow(ex)
            end
        else
            rethrow(ex)
        end
    end
    error("huh? Shouldn't get here")
end


""" 

Solve f(x) = 0 if we have a bracket [a,b] and starting point a < c < b

Use Algorithm 4.2 of Alefeld, Potra, Shi unless it dithers, in which case the root is found by bisection.

"""
function have_bracket(f, a, b, c=(0.5)*(a+b);
                      xtol::Real=zero(c), xtolrel::Real=zero(c), maxeval=10, verbose=false)
    ## try a42a unless it fails, then go to failsafe
    verbose && println("Have a bracket, switching to bracketing method")
    f(a) == 0.0 && return(a)
    f(b) == 0.0 && return(b)
    if f(a) * f(b) > 0
        throw(ConvergenceFailed("Interval [$a, $b] is not a bracket"))
    end
    a,b = sort([a,b])
    if !(a < c < b)
        throw(ConvergenceFailed("Value $c is not in interval ($a, $b)"))
    end

    try
        x0 = a42a(f, a, b, c, xtol=xtol, maxeval=maxeval, verbose=verbose)
    catch e
        if isa(e, StateConverged)
            return e.x0
        else
            verbose && println("Dithering, switching to bisection method")
            return find_zero(f, a, b, verbose=verbose)
        end
    end
end
            

## assume f(b) > 0
## steffensen of secant line method, but we don't assume a bracket
## we "bend" large trajectories 
function secant_method_no_bracket(f, a, b;
                                  ftol=10*eps(one(float(a))), xtol::Real=10*eps(one(float(a))),
                                  xtolrel::Real=10*eps(one(float(a))),
                                  maxeval::Int=100, verbose::Bool=false)
    
    alpha, beta = a, b ## beta keeps smallest values
    falpha, fbeta = f(alpha), f(beta)

    ## checks -- have a bracket
    if falpha * fbeta < 0
        a = a42(f, alpha, beta, xtol=xtol, maxeval=maxeval, verbose=verbose)
        throw(StateConverged(a))
    end
    
    ## keep fbeta smallest
    if falpha < fbeta
        alpha, beta, falpha, fbeta = beta, alpha, fbeta, falpha
    end

    gamma = Inf

    ctr, quad_ctr = 0, 0

    while fbeta > 0
        if ctr > maxeval
            if abs(f(beta)) < ftol
                throw(StateConverged(beta))
            else
                ## try a higher order method? May get us there
                out = kss5(f, beta, ftol=float(ftol), xtol=float(xtol), reltol=float(xtolrel))
                xn = out.x[end]
                if out.state == :converged && !isnan(xn) && !isinf(xn)
                    throw(StateConverged(xn))
                else
                    throw(Dithering(beta))
                end
            end
        else
            (ctr = ctr + 1)
        end
        ## use steffensen or secant line depending on fbeta's size
        if abs(fbeta) < abs(beta - alpha)
            stepsize = fbeta * fbeta / (f(beta + fbeta) - fbeta)  # steffensen
        else
            stepsize =  fbeta * (beta - alpha)  / (fbeta - falpha) # secant line
        end
        gamma = beta - stepsize

        ## bend in some cases. XXX What is best bending scheme? XXX
        if abs(gamma-beta) >= 5e2 * abs(beta - alpha)
            stepsize = cbrt( (beta - alpha) / (fbeta - falpha) * fbeta)
#            gamma = beta + (1/4) * (beta - alpha) / (fbeta - falpha) * fbeta
            gamma = beta - stepsize
        end

        fgamma = f(gamma)

        
        # try
        #     fgamma = f(gamma) ## what to do with domain error?
        # catch e
        #     if isa(e, DomainError)
        #         error("domain error -- fix me")
        #     end
        # end
        
        if abs(fgamma) <= ftol
            throw(StateConverged(gamma))
        end
        if fgamma < 0
            x0 = have_bracket(f, gamma, beta, maxeval=20, verbose=verbose)
            throw(StateConverged(x0))
        end
        
        ## increasing value, try quadratic
        if fgamma >= fbeta
            quad_ctr = quad_ctr + 1
            if quad_ctr > 5*10
                abs(fbeta) <  ftol ? throw(StateConverged(beta)) : throw(Dithering(beta))
            end
            ## secant line stopped decreasing
            if alpha == gamma || gamma == beta
                throw(PossibleExtremaReached(gamma)) ## poorly titled error
            end

            ## quadratic interpolation to get gamma
            denom = (beta - alpha) * (fbeta - fgamma)  - (beta - fgamma) * (fbeta - falpha)
            if denom == 0.0
                throw(PossibleExtremaReached(alpha))
            end
            gamma = beta -  ((beta - alpha)^2 * (fbeta - fgamma) - (beta - gamma)^2 * (fbeta - falpha))/denom/2

            fgamma = f(gamma)

            if fgamma == 0
                throw(StateConverged(gamma))
            end
            if fgamma < 0
                x0 = have_bracket(f, gamma, beta, maxeval=20, verbose=verbose)
                throw(StateConverged(x0))
            end
        else
            quad_ctr = 0        # reset
        end

        if fgamma > fbeta
            alpha, beta, falpha, fbeta = gamma, beta, fgamma, fbeta
        elseif fgamma < fbeta
            alpha, beta, falpha, fbeta = beta, gamma, fbeta, fgamma
        else
            throw(PossibleExtremaReached(beta))
        end
        verbose && println("x_$ctr: $beta\tf(x_$ctr)=$(f(beta))")
    end
    if abs(f(beta)) < ftol  
        throw(StateConverged(beta))
    else
        throw(ConvergenceFailed("More than $maxeval iterations before convergence"))
    end
end
  
