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
    h = eps(one(real(x1)))^(1//3)
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



# Should a Steffensen step or secant step be taken? If |fx| small enough, take a
# Steffensen step. A steffensen step uses f(x + fx) - fx = f'(x)*fx + f''(ξ) /2 * fx^2
# We keep the error small here as long as |x|f''(x) not huge.
@inline function do_steff_step(x::T, fx::S) where {T, S}
    abs(fx) <=  max(oneunit(S), oneunit(S)*abs(x)/oneunit(T)) * one(T) * 1//1000
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

# return vertex of parabola through (a,fa),(b,fb),(c,fc)
# first time trhough, we have picture of a > b > c; |fa|, |fc| > |fb|, all same sign
function quad_vertex(c,fc,b,fb,a,fa)
    fba = (fb-fa)/(b-a)
    fbc = (fb-fc)/(b-c)

    1/2 * ((a+b) - fba/(fbc - fba)*(c-a))

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
