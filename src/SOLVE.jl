## SOLVE button from the HP-34C
## follows roughly algorithm described http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf
## though some modifications were made.
##
## Goal is to return a value `x` with either:
## * `f(x) == 0.0` or 
## * `f(prevfloat(x)) * f(nextfloat(x)) < 0`.
## if a bracket is found that can be done. 
##
## Otherwise, values with `abs(f(x)) < eps()^(1/5)` are given status of answers
## and others are marked as a possible extrema.
##
## If the algorithm gets stuck, a Dithering error is thrown.
##
## While this is slower than fzero with order 2, 5, 8, or 16 it is 
## more robust to initial condition
function SOLVE(f::Function, x0::Real)
    x0 =  float(x0)
    fx0 = f(x0)
    if f(x0) == 0.0 return x0 end
    try
        ## get initial guess. Use Steffensen guess if reasonable size
        step = max(1/100, min(abs(fx0), abs(x0/100)))
        if fx0 < 0
            g(x) = -f(x)
            secant_method_no_bracket(g, x0 - step, x0)
        else
            secant_method_no_bracket(f, x0 - step, x0)
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
            if abs(f(guess)) < eps(zero(x0))^(1/5) 
                return guess
            else
                throw(ex)
            end
        else
            throw(ex)
        end
    end
    error("huh? Shouldn't get here")
end


type PossibleExtremaReached
    x0
end

type Dithering
    x0
end

## Use steffensen or secant method within a bracket unless it dithers
## in which case it uses bisection method
## could use brent's method, might be faster
## finds value x with f(x) == 0.0 of f(prevfloat(x)) * f(nextfloat(x)) <= 0
function secant_method_bracket(f::Function, a, b)
    alpha, beta, gamma = min(a, b), max(a,b), (a + b)/2
    falpha, fbeta = f(alpha), f(beta)
    m, M = alpha, beta

    @assert falpha * fbeta < 0
    
    ctr = 0
    while nextfloat(m) < M
        ## give secant 10 tries, then go failsafe
        ctr > 10 ? throw(StateConverged(fzero(f, [m, M]))) : (ctr = ctr + 1)
        if abs(fbeta) < abs(beta - alpha)
            step = fbeta*fbeta/(f(beta + fbeta) - fbeta)
        else
            step = (beta - alpha) * fbeta / (fbeta - falpha)
        end
        ## gamma is secant line guess or steffensen
        gamma = beta - step
        if gamma <= m || gamma >= M
            gamma = (alpha + beta) / 2 ## "bend via bracket"
        end
        if (gamma == alpha || gamma == beta)
            throw(StateConverged(gamma) )
        end

        fgamma = f(gamma)
        if fgamma == 0.0
            throw(StateConverged(gamma))
        end

        if fgamma * f(M) < 0
            m, M = gamma, M
        else
            m, M = m, gamma
        end
        alpha, beta, falpha, fbeta = beta, gamma, fbeta, fgamma

    end
    throw(StateConverged(gamma))
end

## assume f(b) > 0
## steffensen of secant line method, but we don't assume a bracket
## we "bend" large trajectories 
function secant_method_no_bracket(f::Function, a, b)
    alpha, beta = a, b ## beta keeps smallest values
    falpha, fbeta = f(alpha), f(beta)

    ## checks -- have a bracket
    if falpha * fbeta < 0
        secant_method_bracket(f, alpha, beta)
    end
    
    ## keep fbeta smallest
    if falpha < fbeta
        alpha, beta, falpha, fbeta = beta, alpha, fbeta, falpha
    end

    gamma, delta = Inf, (alpha + beta)/2

    ctr, quad_ctr = 0, 0
    while fbeta > 0
        if ctr > 1000
            if abs(f(beta)) < eps(zero(b))^(1/5) 
                throw(StateConverged(beta))
            else
                throw(Dithering(beta))
            end
        else
            (ctr = ctr + 1)
        end
        ## use steffensen or secant line depending on fbeta's size
        if abs(fbeta) < abs(beta - alpha)
            step = fbeta / (f(beta + fbeta) - fbeta) * fbeta # steffensen
        else
            step =  (beta - alpha) * fbeta / (fbeta - falpha) # secant line
        end
        gamma = beta - step

        ## bend in some cases. What is best bending scheme?
        if abs(gamma-beta) >= 500 * abs(beta - alpha)
#            gamma = beta + (1/4) * (beta - alpha) / (fbeta - falpha) * fbeta
            gamma = beta - cbrt( (beta - alpha) / (fbeta - falpha) * fbeta)
        end

        fgamma = f(gamma)
        # try
        #     fgamma = f(gamma) ## what to do with domain error?
        # catch e
        #     if isa(e, DomainError)
        #         error("domain error -- fix me")
        #     end
        # end
        
        if fgamma == 0
            throw(StateConverged(gamma))
        end
        if fgamma < 0
            secant_method_bracket(f, gamma, beta)
        end
        
        ## increasing value, try quadratic
        if fgamma >= fbeta
            quad_ctr = quad_ctr + 1
            if quad_ctr > 5*10
                abs(fbeta) <  eps(zero(b))^(1/5) ? throw(StateConverged(beta)) : throw(Dithering(beta))
            end
            ## secant line stopped decreasing
            if alpha == gamma || gamma == beta
                throw(PossibleExtremaReached(gamma)) ## poorly titled error
            end

            ## quadratic interpolation to get gamma
            den = (beta - alpha) * (fbeta - fgamma)  - (beta - fgamma) * (fbeta - falpha)
            if den == 0.0
                throw(PossibleExtremaReached(alpha))
            end
            gamma = beta -  ((beta - alpha)^2 * (fbeta - fgamma) - (beta - gamma)^2 * (fbeta - falpha))/den/2

            fgamma = f(gamma)
            if fgamma == 0
                throw(StateConverged(gamma))
            end
            if fgamma < 0
                secant_method_bracket(f, gamma, beta)
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
    end
    throw(StateConverged(beta))
end
  
