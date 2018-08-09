##################################################

# type to throw on succesful convergence
mutable struct StateConverged
    x0::Number
end

# type to throw on failure
mutable struct ConvergenceFailed
    reason::AbstractString
end

##################################################
## Helpers for the various methods

_unitless(x) = x / oneunit(x)

## issue with approx derivative
isissue(x) = iszero(x) || isnan(x) || isinf(x)

# of (a,fa), (b,fb) choose pair where |f| is smallest
@inline choose_smallest(a, b, fa, fb) = abs(fa) < abs(fb) ? (a,fa) : (b,fb)
@inline sort_smallest(a, b, fa, fb) = abs(fa) < abs(fb) ? (a, b, fa, fb) : (b, a, fb, fa)


## find a default secant step
function _default_secant_step(x1)  
    h = eps(one(real(x1)))^(1/3)
    dx = h*oneunit(x1) + abs(x1)*h^2 # adjust for if eps(x1) > h
    x0 = x1 + dx
    x0
end


"""
heuristic to get a decent first step with Steffensen steps
"""
function steff_step(x, fx)

    xbar, fxbar = real(x/oneunit(x)), fx/oneunit(fx)
    thresh =  max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6
    
    out = abs(fxbar) <= thresh ? fxbar  : sign(fx) * thresh 
    out * oneunit(x)
    
end

function guarded_secant_step(alpha, beta, falpha, fbeta)

    fp = (fbeta - falpha) /  (beta - alpha)
    Δ = fbeta / fp
    ## odd, we get allocations if we define Delta, then beta - Delta
    ## Δ = beta - fbeta * (beta - alpha) / (fbeta - falpha)

    if isissue(Δ)
        Δ = oneunit(alpha)/1000
    elseif abs(Δ) >= 100 * abs(alpha - beta) # guard runaway
        Δ = sign(Δ) * 100 * min(oneunit(alpha), abs(alpha - beta))
    end

    if isissue(Δ)
        return (alpha + (beta - alpha)*(0.5), true) # midpoint
    else
        return (beta - Δ, false)
    end
end


# for the 3 points, find parabola. Then return vertex of the closest 0
# to the mean of x1, x2, x3
function quad_vertex(x1,fx1,x2,fx2,x3,fx3)
    vertex = -(-fx1*(x2^2 - x3^2) + fx2*(x1^2 - x3^2) - fx3*(x1^2 - x2^2))/(2*(fx1*(x2 - x3) - fx2*(x1 - x3) + fx3*(x1 - x2)))
    discr = -4*(fx1*(x2 - x3) - fx2*(x1 - x3) + fx3*(x1 - x2))*(fx1*x2*x3*(x2 - x3) - fx2*x1*x3*(x1 - x3) + fx3*x1*x2*(x1 - x2))/(x1^2*x2 - x1^2*x3 - x1*x2^2 + x1*x3^2 + x2^2*x3 - x2*x3^2)^2 + (-fx1*(x2^2 - x3^2) + fx2*(x1^2 - x3^2) - fx3*(x1^2 - x2^2))^2/(x1^2*x2 - x1^2*x3 - x1*x2^2 + x1*x3^2 + x2^2*x3 - x2*x3^2)^2

    if discr > zero(discr)
        b = (-fx1*(x2^2 - x3^2) + fx2*(x1^2 - x3^2) - fx3*(x1^2 - x2^2))/(x1^2*x2 - x1^2*x3 - x1*x2^2 + x1*x3^2 + x2^2*x3 - x2*x3^2)
        a = (fx1*(x2 - x3) - fx2*(x1 - x3) + fx3*(x1 - x2))/(x1^2*x2 - x1^2*x3 - x1*x2^2 + x1*x3^2 + x2^2*x3 - x2*x3^2)
        gamma1, gamma2 = (-b + sqrt(discr))/2a, (-b - sqrt(discr))/2a

        xbar = (x1+x2+x3)/3
        d1,d2,d3 = abs(gamma1 - xbar), abs(gamma2-xbar), abs(vertex - xbar)
        # return closest
        d1 < min(d2,d3) && return gamma1
        d2 < min(d1,d3) && return gamma2
    end
    return vertex
end



## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    num, den = fb - fa, b - a
    iszero(num) && iszero(den) && return Inf, true
    out = num / den
    out, isissue(out)
end

## use f[y,z] - f[x,y] + f[x,z] to approximate
function _fbracket_diff(a,b,c, fa, fb, fc)
    x1, issue = _fbracket(b, c, fb,  fc)
    issue && return x1, issue
    x2, issue = _fbracket(a, b, fa,  fb)
    issue && return x2, issue    
    x3, issue = _fbracket(a, c, fa,  fc)
    issue && return x3, issue
    
    out = x1 - x2 + x3
    out, isissue(out)
end


## use f[a,b] * f[a,c] / f[b,c]
function _fbracket_ratio(a, b, c, fa, fb, fc)
    x1, _ = _fbracket(a, b, fa, fb)
    x2, _ = _fbracket(a, c, fa, fc)
    x3, _ = _fbracket(b, c, fb, fc)
    out = (x2 * x3) / x3
    out, isissue(out)
end

## from https://core.ac.uk/download/pdf/82282744.pdf
## signum based iteration allows one to guarantee a valid starting point
## if N is big enough. (THough this depends on algorithm, a, b and function)
## N here would need to be tuned. But, as is, this may prove useful anyways.
function identify_starting_point(f, a, b, N)
    pts = range(a, stop=b, length=N+1)
    fxs = f.(pts)
    sfxs = sign.(f.(pts))
    identify_starting_point(a,b,sfxs)
end

function identify_starting_point(a, b, sfxs)
    N = length(sfxs) - 1
    p0 = a + (b-a)/2
    p1 = p0 + (b-a)/(2N) * sfxs[1] * sum(s for s in sfxs[2:end-1])
    p1
end
