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
struct Order0 <: AbstractSecant end

""" find_zero(y, 1.8s, order) 
    Order1()

The `Order1()` method is an alias for `Secant`. It specifies the
[secant method](https://en.wikipedia.org/wiki/Secant_method).
This method keeps two values in its state, `x_n` and `x_n1`. The
updated point is the intersection point of x axis with the secant line
formed from the two points. The secant method uses 1 function
evaluation per step and has order `(1+sqrt(5))/2`.

"""
struct Secant <: AbstractSecant end
const Order1 = Secant


function init_state(method::AbstractSecant, fs, x::Union{Tuple, Vector})
    x0, x1 = promote(float(x[1]), float(x[2]))
    fx0, fx1 = fs(x0), fs(x1)        
    UnivariateZeroState(x1, x0,
                        missing, #oneunit(x1) * (0*x1)/(0*x1),
                        
                        fx1, fx0, 0, 2,
                        false, false, false, false, "")
end

function init_state(method::AbstractSecant, fs, x::Number)

    # need an initial x0,x1 if two not specified
    x1 = float(x)
    fx1 = fs(x1)

    h = eps(one(real(x1)))^(1/3)
    dx = h*oneunit(x1) + abs(x1)*h^2 # adjust for if eps(x1) > h
    x0 = x1 + dx
    fx0 = fs(x0)
    
    UnivariateZeroState(x1, x0,
                        missing, #oneunit(x1) * (0*x1)/(0*x1),
                        fx1, fx0, 0, 2,
                        false, false, false, false, "")

end

function init_state!(state::UnivariateZeroState{T, S}, ::AbstractSecant, f, x::Union{Tuple, Vector}) where {T, S}
    x0,x1 = promote(float.(x))
    fx0, fx1 = promote(f(x0), f(x1))
    init_state!(state, x0, x1, missing, fx0, fx1)
end

function init_state!(state::UnivariateZeroState{T, S}, ::AbstractSecant, fs, x::Number) where {T, S}
    x1 = float(x)
    h = eps(one(real(x1)))^(1/3)
    dx = h*oneunit(x1) + abs(x1)*h^2 # adjust for if eps(x1) > h
    x0 = x1 + dx
    fx0, fx1 = promote(fs(x0), fs(x1))
    
    init_state!(state, x0, x1, missing, fx0, fx1)
end
##################################################


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

## in Order0, we run bisection if a bracketing interval is found
## this is meant to be as speedy as possible
function _run_bisection(fs, options, state)
    ## could do hybrid method here, where we update state and options, but
    ## that proves slower, as bisection64 is more efficient.
    xn0, xn1 = state.xn0, state.xn1
    state.xn1 = bisection(fs, xn0, xn1)
    state.x_converged = true
    options.strict = true # prevent check on f(xn)
    state.message *= "Used bisection for last step, evaluation counts are not accurate. "
end

function update_state(method::Order0, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    alpha, beta = o.xn0, o.xn1
    falpha, fbeta = o.fxn0, o.fxn1

    if sign(falpha) * sign(fbeta) < 0.0
        _run_bisection(fs, options, o)
        return nothing
    end

    gamma::T, issue = guarded_secant_step(alpha, beta, falpha, fbeta)

    fgamma::S = fs(gamma); incfn(o)

    if sign(fgamma)*sign(fbeta) < 0.0
        o.xn0, o.xn1 = gamma, beta
        o.fxn0, o.fxn1 = fgamma, fbeta
        _run_bisection(fs, options, o)
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
        gamma = quad_vertex(alpha, falpha, beta, fbeta, gamma, fgamma)
#        denom = (beta - alpha) * (fbeta - fgamma)  - (beta - fgamma) * (fbeta - falpha)
#        if isissue(denom)
        if isissue(gamma)            
            o.stopped
            o.message = "dithering, algorithm failed to improve using quadratic steps"
            return nothing
        end
        #        gamma = beta -  ((beta - alpha)^2 * (fbeta - fgamma) - (beta - gamma)^2 * (fbeta - falpha))/denom/2



        fgamma = fs(gamma); incfn(o)
        incfn(o)

        if abs(fgamma) < abs(fbeta)
            o.xn0, o.xn1 = beta, gamma
            o.fxn0, o.fxn1 = fbeta, fgamma
            return nothing
        end

        theta::T, issue = guarded_secant_step(beta, gamma, fbeta, fgamma)
        
        ftheta::S = fs(theta); incfn(o)

        if sign(ftheta) * sign(fbeta) < 0
            o.xn0, o.xn1 = beta, theta
            o.fxn0, o.fxn1 = fbeta, ftheta

            _run_bisection(fs, options, o)
            return nothing
        end
    end

    # failed to improve
    o.stopped = true
    o.message = "failure to improve"
    return nothing
end

function show_tracks(l::Tracks, M::Order0)
    println("Tracks not recorded for Order0()")
    println("")
end

##################################################

## Secant
## https://en.wikipedia.org/wiki/Secant_method

function update_state(method::Secant, fs, o::UnivariateZeroState{T,S}, options)  where {T, S}

    if (o.fxn0 == o.fxn1) || (o.xn0 == o.xn1) 
         o.stopped = true
         o.message = "Derivative approximation had issues"
         return
     end

    dx = o.fxn1 * (o.xn1 - o.xn0) / (o.fxn1 - o.fxn0)
    o.xn0, o.xn1 = o.xn1, o.xn1 - dx
    o.fxn0, o.fxn1 = o.fxn1, fs(o.xn1)
    incfn(o)
    
    nothing

end



##################################################

### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description
struct Steffensen <: AbstractUnivariateZeroMethod
end

"""
    Order2()

The quadratically converging
[Steffensen](https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description)
method is used for the derivative-free `Order2()` algorithm. Unlike
the quadratically converging Newton's method, no derivative is
necessary, though like Newton's method, two function calls per step
are. This algorithm is more sensitive than Newton's method to poor
initial guesses.

"""    
const Order2 = Steffensen

function update_state(method::Steffensen, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    
    wn::T = o.xn1 + steff_step(o.xn1, o.fxn1)

    fwn::S = fs(wn)
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

Implements an order 5 algorithm from *A New Fifth Order Derivative
Free Newton-Type Method for Solving Nonlinear Equations* by Manoj
Kumar, Akhilesh Kumar Singh, and Akanksha, Appl. Math. Inf. Sci. 9,
No. 3, 1507-1513 (2015), DOI: 10.12785/amis/090346. Four function
calls per step are needed.

"""    
struct Order5 <: AbstractUnivariateZeroMethod end


## If we have a derivative, we have this
function update_state(method::Order5, fs::Union{FirstDerivative,SecondDerivative},
                      o::UnivariateZeroState{T,S}, options)  where {T, S}


    xn, fxn = o.xn1, o.fxn1

    fpxn::S = fs(xn, 1)
    incfn(o)

    if isissue(fpxn)
        o.stopped  = true
        return
    end

    yn::T = xn - fxn / fpxn
    fyn::S, fpyn::S = fs(yn), fs(yn, 1)
    incfn(o, 2)

    if isissue(fpyn)
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.stopped  = true
        return
    end


    zn::T = xn  - (fxn + fyn) / fpxn
    fzn::S = fs(zn)
    incfn(o, 1)

    xn1 = zn - fzn / fpyn
    fxn1 = fs(xn1)
    incfn(o, 1)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


function update_state(method::Order5, fs, o::UnivariateZeroState{T,S}, options) where {T, S}

    xn = o.xn1
    fxn = o.fxn1

    wn::T = o.xn1 + steff_step(o.xn1, o.fxn1)

    fwn::S = fs(wn)
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
    fyn::S = fs(yn)
    incfn(o)


    zn::T = xn - (fxn + fyn) / fp
    fzn::S = fs(zn)
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


## cf also: https://doi.org/10.1515/tmj-2017-0049
"""
    Order8()

Implements an eighth-order algorithm from *New Eighth-Order
Derivative-Free Methods for Solving Nonlinear Equations* by Rajinder
Thukral, International Journal of Mathematics and Mathematical
Sciences Volume 2012 (2012), Article ID 493456, 12 pages DOI:
10.1155/2012/493456. Four function calls per step are required.
    
"""
struct Order8 <: AbstractUnivariateZeroMethod
end

function update_state(method::Order8, fs, o::UnivariateZeroState{T,S}, options) where {T, S}

    xn = o.xn1
    fxn = o.fxn1

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
struct Order16 <: AbstractUnivariateZeroMethod
end

function update_state(method::Order16, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    xn = o.xn1
    fxn = o.fxn1

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


