# Many derivative free methods of different orders

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

function init_state(method::AbstractSecant, fs, x)

    if isa(x, Vector) || isa(x, Tuple)
        x0, x1 = x[1], x[2]
        fx0, fx1 = fs(x0), fs(x1)        
    else
        # need an initial x0,x1 if two not specified
        x0 = x
        fx0 = fs(x0)        
        stepsize = max(1/100, min(abs(fx0/oneunit(fx0)), abs(x0/oneunit(x0)/100)))
        x1 = x0 + stepsize*oneunit(x0)
        x0, x1, fx0, fx1  = x1, x0, fs(x1), fx0 # switch        
    end

    state = UnivariateZeroState( promote(x1, x0)...,
                                 promote(fx1, fx0)...,
                                 0, 2,
                                 false, false, false, false,
                                 "")
    state

end

##################################################

## in Order0, we run bisection if a bracketing interval is found
## with verbose turned off
function _run_bisection(fs, o, options)
    verbose = options.verbose; options.verbose=false # turn off verbose
    find_zero(Bisection(), fs, options, o)
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
function update_state(method::Order0, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions) where {T}

    alpha, beta = o.xn0, o.xn1
    falpha, fbeta = o.fxn0, o.fxn1

    incsteps(o)

    if sign(falpha) * sign(fbeta) < 0.0
        _run_bisection(fs, o, options)
        return nothing
    end

    gamma, issue = guarded_secant_step(alpha, beta, falpha, fbeta)

    fgamma = fs(gamma); incfn(o)

    if sign(fgamma)*sign(fbeta) < 0.0
        o.xn0, o.xn1 = gamma, beta
        o.fxn0, o.fxn1 = fgamma, fbeta
        _run_bisection(fs, o, options)
        return nothing
    end

    if abs(fgamma) <= abs(fbeta)
        o.xn0, o.xn1 = beta, gamma
        o.fxn0, o.fxn1 = fbeta, fgamma
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

        # quadratic_step. Put new gamma at vertex of parabola through alpha, beta, (old) gamma
        denom = (beta - alpha) * (fbeta - fgamma)  - (beta - fgamma) * (fbeta - falpha)
        if isissue(denom)
            o.stopped
            o.message = "dithering, algorithm failed to improve using quadratic steps"
            return nothing
        end
        gamma = beta -  ((beta - alpha)^2 * (fbeta - fgamma) - (beta - gamma)^2 * (fbeta - falpha))/denom/2


        fgamma = fs(gamma); incfn(o)
        incfn(o)

        if abs(fgamma) < abs(fbeta)
            o.xn0, o.xn1 = beta, gamma
            o.fxn0, o.fxn1 = fbeta, fgamma
            return nothing
        end

        theta, issue = guarded_secant_step(beta, gamma, fbeta, fgamma)
        
        ftheta = fs(theta); incfn(o)

        if sign(ftheta) * sign(fbeta) < 0
            o.xn0, o.xn1 = beta, theta
            o.fxn0, o.fxn1 = fbeta, ftheta

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
function update_state(method::Secant, fs, o, options) 

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
    o.fxn1 = fs(o.xn1)    
    incfn(o)

    nothing
end

"""

    secant_method(f, x0, x1; [kwargs...])
    
Solve for zero of `f(x) = 0` using the secant method.

Not exported.  Use `find_zero` with `Order1()`.    
"""
secant_method(f, x0::Number, x1::Number; kwargs...) = find_zero(f, (x0, x1), Order1(); kwargs...)


##################################################

### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description
mutable struct Steffensen <: AbstractUnivariateZeroMethod
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

function update_state(method::Steffensen, fs, o, options) 

    incsteps(o)

    wn = o.xn1 + steff_step(o.xn1, o.fxn1)
    fwn = fs(wn)
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
mutable struct Order5 <: AbstractUnivariateZeroMethod end

## If we have a derivative, we have this
function update_state(method::Order5, fs::Union{FirstDerivative,SecondDerivative}, o, options) 


    xn, fxn = o.xn1, o.fxn1

    incsteps(o)

    fpxn = fs(xn, 1)
    incfn(o)

    if isissue(fpxn)
        o.stopped  = true
        return
    end

    yn = xn - fxn / fpxn
    fyn, fpyn = fs(yn), fs(yn, 1)
    incfn(o, 2)

    if isissue(fpyn)
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped  = true
        return
    end


    zn = xn  - (fxn + fyn) / fpxn
    fzn = fs(zn)
    incfn(o, 1)

    xn1 = zn - fzn / fpyn
    fxn1 = fs(xn1)
    incfn(o, 1)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


function update_state(method::Order5, fs, o, options)

    xn = o.xn1
    fxn = o.fxn1

    incsteps(o)

    wn = o.xn1 + steff_step(o.xn1, o.fxn1)

    fwn = fs(wn)
    incfn(o)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)
    if issue
        o.xn0, o.xn1 = o.xn1, wn
        o.fxn0, o.fxn1 = o.fxn1, fwn
        o.message = "Issue with divided difference f[xn, wn]"
        o.stopped  = true
        return
    end

    yn = o.xn1 - o.fxn1 / fp
    fyn = fs(yn)
    incfn(o)


    zn = xn - (fxn + fyn) / fp
    fzn = fs(zn)
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

##################################################



"""
    Order8()

Implements an algorithm from 
*New Eighth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by Rajinder Thukral,
International Journal of Mathematics and Mathematical Sciences
Volume 2012 (2012), Article ID 493456, 12 pages. Four function calls per step are required.
"""
mutable struct Order8 <: AbstractUnivariateZeroMethod
end

function update_state(method::Order8, fs, o, options) 

    xn = o.xn1
    fxn = o.fxn1
    incsteps(o)
    S = eltype(fxn)

    wn = xn + steff_step(xn, fxn)
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

    yn = xn - fxn / fp
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

    w = 1 / (1 - fzn/fwn)

    xi = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))

    xn1 = zn - w * xi * fzn / fp
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

Five function calls per step are required. Though rapidly converging,
this method generally isn't faster (fewer function calls/steps) over
other methods when using `Float64` values, but may be useful for
solving over `BigFloat`.

"""
mutable struct Order16 <: AbstractUnivariateZeroMethod
end

function update_state(method::Order16, fs, o, options) 
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn = xn + steff_step(xn, fxn)
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

    yn = xn - fxn / fp
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



    zn = yn - fyn / fp
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
    fan::S = fs(an)
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

    fp1, issue = _fbracket(xn,yn, fxn, fyn)
    
    sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
    u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/(fp1/oneunit(fp1))


    xn1 = an - sigma * fan / fp
    fxn1::S = fs(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


