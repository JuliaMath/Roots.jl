# Many derivative free methods of different orders

##################################################
## Helpers for the various methods
## issue with approx derivative
isissue(x) = (x == 0.0) || isnan(x) || isinf(x)


"""
heuristic to get a decent first step with Steffensen steps
"""
function steff_step(x::T, fx) where {T}
    thresh =  max(1, norm(x)) * sqrt(eps(T)) # max(1, sqrt(abs(x/fx))) * 1e-6
    norm(fx) <= thresh ? fx : sign(fx) * thresh
end


function guarded_secant_step(alpha::T, beta::T, falpha, fbeta) where {T <: AbstractFloat}


    fp = (fbeta - falpha) /  (beta - alpha)
    Δ::T = fbeta / fp
    ## odd, we get allocations if we define Delta, then beta - Delta
    ## Δ = beta - fbeta * (beta - alpha) / (fbeta - falpha)
    
    if isissue(Δ)
        Δ = one(T)/1000
    elseif norm(Δ) >= 100 * norm(alpha - beta) # guard runaway
        Δ = sign(Δ) * 100 * min(one(T), norm(alpha - beta))
    end


    if isissue(Δ)
        return (alpha + (beta - alpha)*(0.5), true) # midpoint
    else
        return (beta - Δ, false)
    end


end


## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    num, den = fb - fa, b - a
    num == 0 && den == 0 && return Inf, true
    out = num / den
    out, isissue(out)
end

## use f[y,z] - f[x,y] + f[x,z] to approximate
function _fbracket_diff(a,b,c, fa, fb, fc)
    x1, _ = _fbracket(b, c, fb,  fc)
    x2, _ = _fbracket(a, b, fa,  fb)
    x3, _ = _fbracket(a, c, fa,  fc)
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


##################################################

## Order0 and Secant are related
"""
    Order0()


The `Order0` method is engineered to be a more robust, though possibly
slower, alternative to to the other derivative-free root-finding
methods. The implementation roughly follows the algorithm described in
*Personal Calculator Has Key to Solve Any Equation f(x) = 0*, the
SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
The basic idea is to use a secant step. If along the way a bracket is
found, switch to bisection, using `Bisection` if possible, else
`A42`.  If the secant step fails to decrease the function value, a
quadratic step is used up to 3 times.
    
"""
mutable struct Order0 <: AbstractSecant end

"""
    Order1()

The `Order1()` method is an alias for `Secant` and performs a secant
step. This method keeps two values in its state, `x_n` and `x_n1`. The
updated point is the intersection point of x axis with the secant line
formed from the two points. The secant method uses 1 function
evaluation and has order `(1+sqrt(5))/2`.

"""
mutable struct Secant <: AbstractSecant end
const Order1 = Secant

function init_state(method::AbstractSecant, fs, x::Union{T, Vector{T}}, bracket) where {T <: AbstractFloat}

    if isa(x, Vector)
        x0, x1 = x[1:2]
        #        fx0, fx1 = fs.f(x0), fs.f(x1)
        fx0, fx1 = fs(x0), fs(x1)        
    else
        x0 = float(x)
        #        fx0 = fs.f(x0)
        fx0 = fs(x0)        
        stepsize = max(1/100, min(abs(fx0), abs(x0/100)))
        x1 = x0 + stepsize
#        x0, x1, fx0, fx1  = x1, x0, fs.f(x1), fx0 # switch
        x0, x1, fx0, fx1  = x1, x0, fs(x1), fx0 # switch        
    end

    state = UnivariateZeroStateBase(
                                    promote(float(x1), float(x0))...,
                                    promote(fx1, fx0)...,
                                    bracket,
                                    0, 2,
                                    false, false, false, false,
                                    "")
    state
end

##################################################

## in Order0, we run bisection if a bracketing interval is found
function _run_bisection(fs, o, options)
    verbose = options.verbose; options.verbose=false # turn off verbose
    find_zero(Bisection(), fs, o, options)
    options.verbose = verbose
    o.message = "Used bisection for last step"
end


## order 0
# goal: more robust to initial guess than higher order methods
# follows roughly algorithm described http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf, the SOLVE button from the HP-34C
# though some modifications were made.
# * use secant step
# * if along the way a bracket is found, switch to bisection. (We use float64 bisection not a42 if available)
# * if secant step fails to decrease, we use quadratic step up to 3 times
#
# Goal is to return a value `x` with either:
# * `f(x) == 0.0` or
# * `f(prevfloat(x)) * f(nextfloat(x)) < 0`.
# if a bracket is found that can be done, otherwise secant step is used
function update_state(method::Order0, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions{T}) where {T}

    f = fs
    α = o.xn0
    β = o.xn1
    fα =  o.fxn0
    fβ = o.fxn1

    S = eltype(fα)

    incsteps(o)

    if sign(fα) * sign(fβ) < 0.0
        _run_bisection(fs, o, options)
        return nothing
    end


    gamma, issue = guarded_secant_step(alpha, beta, falpha, fbeta)


    fξ = fs(ξ)
    incfn(o)

    if sign(fξ) * sign(fβ) < 0.0
        o.xn0 = ξ
        o.xn1 =  β
        o.fxn0 = fξ
        o.fxn1 = fβ
        _run_bisection(fs, o, options)
       return nothing
    end
    
    if norm(fξ) <= norm(fβ)
        o.xn0 = β
        o.xn1 = ξ
        o.fxn0 = fβ
        o.fxn1 = fξ
        return nothing
    end
    
    ctr = 0
    while true
        ## quadratic step
        ctr += 1
        if ctr >= 3
           o.stopped = true
           o.message = "dithering, algorithm failed to improve using quadratic steps"
           return nothing
        end

        # quadratic_step. Put new ξ at vertex of parabola through α, β, (old) ξ
        denom = (β - α) * (fβ - fξ)  - (β - fξ) * (fβ - fα)
        if isissue(denom)
            o.stopped = true
            o.message = "dithering, algorithm failed to improve using quadratic steps"
            return nothing
        end
        ξ = tmp::T = β -  ((β - α)^2 * (fβ - fξ) - (β - ξ)^2 * (fβ - fα))/denom/2


        fξ = ftmp::S = f(ξ); incfn(o)
        incfn(o)

        if norm(fξ) < norm(fβ)
            o.xn0, o.xn1 = β, ξ
            o.fxn0, o.fxn1 = fβ, fξ
            return nothing
        end


        theta, issue = guarded_secant_step(beta, gamma, fbeta, fgamma)
        
        ftheta = f(theta); incfn(o)


        if sign(fθ) * sign(fβ) < 0

            o.xn0 = β
            o.xn1 = θ
            o.fxn0 = fβ
            o.fxn1 = fθ

            _run_bisection(fs, o, options)
            return nothing
        end

    end

    # failed to improve
    o.stopped = true
    o.message = "failure to improve"
    return nothing
end


##################################################

## Secant
## https://en.wikipedia.org/wiki/Secant_method
function update_state(method::Secant, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions{T}) where {T <: AbstractFloat}

    incsteps(o)

    fp, issue = _fbracket(o.xn0, o.xn1, o.fxn0, o.fxn1)
    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end

    o.xn0 = o.xn1
    o.fxn0 = o.fxn1

    o.xn1 = o.xn1 -  o.fxn1 / fp
    #    o.fxn1 = fs.f(o.xn1)
    o.fxn1 = fs(o.xn1)    
    incfn(o)

    nothing
end

"""

    secant_method(f, x0, x1; [kwargs...])
    
Solve for zero of `f(x) = 0` using the secant method.

Not exported.  Use `find_zero` with `Order1()`.    
"""
secant_method(f, x0::Real, x1::Real; kwargs...) = find_zero(f, map(float, [x0,x1]), Order1(); kwargs...)


##################################################

### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description
mutable struct Steffensen <: UnivariateZeroMethod
end

"""
    Order2()

The quadratically converging
[Steffensen](https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description)
method is used for the derivative free `Order2()` algorithm. Unlike
the quadratically converging Newton's method, no derivative is
necessary, though like Newton's method, two function calls per step
are. This algorithm is more sensitive than Newton's method to poor
initial guesses.

"""    
const Order2 = Steffensen

function update_state(method::Steffensen, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions{T}) where {T <: AbstractFloat}

    S = eltype(o.fxn1)

    incsteps(o)

    wn = o.xn1 + steff_step(o.xn1, o.fxn1)::T
    fwn = fs(wn)::S
    incfn(o)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)

    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end

    o.xn0 = o.xn1
    o.fxn0 = o.fxn1
    o.xn1 = o.xn1 - o.fxn1 / fp #xn1
    o.fxn1 = fs(o.xn1)
    incfn(o)


    nothing
end

steffenson(f, x0; kwargs...) = find_zero(f, x0, Steffensen(); kwargs...)

##################################################


"""
    Order5()

Implements an algorithm
from *A New Fifth Order Derivative Free Newton-Type Method for Solving Nonlinear Equations*
by Manoj Kumar, Akhilesh Kumar Singh, and Akanksha,
Appl. Math. Inf. Sci. 9, No. 3, 1507-1513 (2015). Four function calls per step are needed.

"""    
mutable struct Order5 <: UnivariateZeroMethod end

function update_state(method::Order5, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}

    xn = o.xn1
    fxn = o.fxn1
    S = eltype(o.fxn1)

    incsteps(o)

    wn::T = o.xn1 + steff_step(o.xn1, o.fxn1)
    fwn = fs(wn)::S
    incfn(o)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)
    if issue
        o.xn0, o.xn1 = o.xn1, wn
        o.fxn0, o.fxn1 = o.fxn1, fwn
        o.message = "Issue with divided difference f[xn, wn]"
        o.stopped  = true
        return
    end

    yn::T = o.xn1 - o.fxn1 / fp
    fyn = fs(yn)::S
    incfn(o)


    zn::T = xn - (fxn + fyn) / fp
    fzn = fs(zn)::S
    incfn(o)

    fp, issue = _fbracket_ratio(yn, o.xn1, wn, fyn, o.fxn1, fwn)
    if issue
        o.xn0, o.xn1 = o.xn1, yn
        o.fxn0, o.fxn1 = o.fxn1, fyn
        o.message = "Issue with f[xn,yn]*f[yn,wn] / f[xn, wn]"
        o.stopped = true
        return
    end

    o.xn0 = o.xn1
    o.fxn0 = o.fxn1
    o.xn1 = zn  - fzn  / fp
    o.fxn1 = fs(o.xn1)
    incfn(o)

    nothing
end

## If we have a derivative
function update_state(method::Order5, fs::FirstDerivative, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}


    xn, fxn = o.xn1, o.fxn1
    S = eltype(fxn)

    incsteps(o)

    fpxn = fs.fp(xn)
    incfn(o)

    if isissue(fpxn)
        o.stopped  = true
        return
    end

    yn = xn - fxn / fpxn
    fyn, fpyn = fs.f(yn), fs.fp(yn)
    incfn(o, 2)

    if isissue(fpyn)
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped  = true
        return
    end


    zn = xn  - (fxn + fyn) / fpxn
    fzn = fs.f(zn)
    incfn(o, 1)

    xn1 = zn - fzn / fpyn
    fxn1 = fs.f(xn1)
    incfn(o, 1)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


##################################################



"""
    Order8()

Implements an algorithm from 
*New Eighth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by Rajinder Thukral,
International Journal of Mathematics and Mathematical Sciences
Volume 2012 (2012), Article ID 493456, 12 pages. Four function calls per step are required.
"""
mutable struct Order8 <: UnivariateZeroMethod
end

function update_state(method::Order8, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs(wn)
    incfn(o)

    if isissue(fwn)
        o.xn0,o.xn1 = xn, wn
        o.fxn0,o.fxn1 = fxn, fwn
        o.stopped = true
        o.message = "issue with Steffensen step fwn"
        return
    end



    fp, issue = _fbracket(xn, wn, fxn, fwn)
    issue && return (xn, true)

    if issue
        o.stopped = true
        o.message = "issue with divided difference f[xn, wn]"
        return
    end

    yn::T = xn - fxn / fp
    fyn::S = fs(yn)
    incfn(o)

    fp, issue = _fbracket(yn, xn, fyn, fxn)
    if issue #fp
        o.xn0,o.xn1 = xn, yn
        o.fxn0,o.fxn1 = fxn, fyn
        o.stopped = true
        o.message = "issue with divided difference f[xn, yn]"
        return
    end


    phi = (1 + fyn / fwn)           # pick one of options
    zn =  yn - phi * fyn / fp
    fzn::S = fs(zn)
    incfn(o)

    fp, issue =  _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    if issue
        o.xn0,o.xn1 = xn, zn
        o.fxn0,o.fxn1 = fxn, fzn
        o.message = "issue with divided difference  f[y,z] - f[x,y] + f[x,z]"
        o.stopped = true
        return
    end

    w::T = 1 / (1 - fzn/fwn)

    xi::T = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))

    xn1::T = zn - w * xi * fzn / fp
    fxn1::S = fs(xn1)
    incfn(o)

    o.xn0,o.xn1 = xn, xn1
    o.fxn0,o.fxn1 = fxn, fxn1

    nothing
end

##################################################

"""
    Order16()

Implement the algorithm from
*New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by R. Thukral, 
American Journal of Computational and Applied Mathematics
p-ISSN: 2165-8935;    e-ISSN: 2165-8943; 2012;  2(3): 112-118
doi: 10.5923/j.ajcam.20120203.08.

Five function calls per step are required. Though rapidly converging, this method generally isn't faster (fewer
function calls/steps) over other methods when using `Float64` values,
but may be useful for solving over `BigFloat`.
"""
mutable struct Order16 <: UnivariateZeroMethod
end

function update_state(method::Order16, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs(wn)
    incfn(o)

    fp, issue = _fbracket(xn, wn, fxn, fwn)

  
    if issue
        o.xn0, o.xn1 = xn, wn
        o.fxn0, o.fxn1 = fxn, fwn
        o.message = "issue with f[xn,wn]"
        o.stopped = true
        return
    end

    yn::T = xn - fxn / fp
    fyn::S = fs(yn)
    incfn(o)

    fp, issue = _fbracket_ratio(yn, xn, wn, fyn, fxn, fwn)
    if issue
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.message = "issue with f[xn,yn]*f[yn,wn]/f[xn,wn]"
        o.stopped = true
        return
    end



    zn::T = yn - fyn / fp
    fzn::S = fs(zn)
    incfn(o)

        
    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn
    eta = 1 / (1 + 2*u3*u4^2) / (1 - u2)
    if issue
        o.xn0, o.xn1 = xn, zn
        o.fxn0, o.fxn1 = fxn, fzn
        o.stopped = true
        o.message = "Approximate derivative failed"
        return
    end

    an = zn - eta * fzn / fp
    fan = fs(an)
    incfn(o)

        
    fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
    if issue
        o.xn0, o.xn1 = xn, an
        o.fxn0, o.fxn1 = fxn, fan
        o.stopped = true
        o.message = "Approximate derivative failed"
        return
    end

    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
    sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
    u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/_fbracket(xn,yn, fxn, fyn)[1]



    

    xn1 = an - sigma * fan / fp
    fxn1 = fs(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


