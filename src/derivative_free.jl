# Many derivative free methods of different orders
#
# TODO: rework Order5 #https://pdfs.semanticscholar.org/ce50/3210d96f653a14b28da96600d5990d2abe97.pdf
# https://content.sciendo.com/view/journals/tmj/10/4/article-p103.xml 7 and 8
# order8: https://www.hindawi.com/journals/ijmms/2012/493456/ref/




## Order0 and 1 are secant type
function init_state(method::AbstractSecant, fs, x::Number)
    x1 = float(x)
    x0 = _default_secant_step(x1)
    init_state(method, fs, (x0, x1))
end


function init_state(method::AbstractSecant, fs, x::Union{Tuple, Vector})
    x0, x1 = promote(float(x[1]), float(x[2]))
    fx0, fx1 = fs(x0), fs(x1)        
    state = UnivariateZeroState(x1, x0, eltype(x1)[],
                                fx1, fx0, eltype(fx1)[],
                                0, 2,
                                false, false, false, false, "")

    state
end



function init_state!(state::UnivariateZeroState{T, S}, method::AbstractSecant, fs, x::Number) where {T, S}
    x1::T = float(x)
    x0::T = _default_secant_step(x1)
    init_state!(state, method, fs, (x0, x1))
end


function init_state!(state::UnivariateZeroState{T, S}, ::AbstractSecant, f, x::Union{Tuple, Vector}) where {T, S}
    x0,x1 = promote(float.(x)...)
    fx0, fx1 = promote(f(x0), f(x1))
    init_state!(state, x1, x0, T[], fx1, fx0, S[])
    state.fnevals = 2
    nothing
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
found, switch to bisection, using `AlefeldPotraShi`.  If the secant
step fails to decrease the function value, a quadratic step is used up
to 3 times.

"""
struct Order0 <: AbstractSecant end

##################################################


## order 0
# goal: more robust to initial guess than higher order methods
# follows roughly algorithm described http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf, the SOLVE button from the HP-34C
# though some modifications were made.
# * use secant step
# * if along the way a bracket is found, switch to bracketing method.
# * if secant step fails to decrease, we use quadratic step up to 3 times
#
# Goal *was* to return a value `x` with either. This can be done if `Bisection` is used
# but for now we opt for a bit speedier using AlefeldPotraShi
# * `f(x) == 0.0` or
# * `f(prevfloat(x)) * f(nextfloat(x)) < 0`.
# if a bracket is found that can be done, otherwise secant step is used
# init_state/init_options try to call Order1 (AbstractSecant) methods with modifications
# we slip in quad_ctr into a value and keep xn1 with smaller norm
function init_state(M::Order0, f, x::Union{Tuple, Vector})
    x0, x1 = promote(float.(x)...)
    fx0, fx1 = promote(f(x0), f(x1))
    
    # we keep xn1, fxn1 smallest
    a,b,fa, fb = sort_smallest(x0, x1, fx0, fx1)


    T, S = eltype(x1), eltype(fx1)
    quad_ctr = one(x1)


    state = UnivariateZeroState(b, a, [one(T)], ## x1, x0, quad_ctr
                                fb, fa, S[], ## fx1, fx0, fc, mflag
                                0, 2,
                                false, false, false, false,
                                "")
    state
end

function init_state!(state::UnivariateZeroState{T,S}, M::Order0, f, x::Union{Tuple,Vector}) where {T,S}
    x0::T, x1::T = x
    fx0::S, fx1::S = promote(f(x0), f(x1))
    a,b,fa, fb = sort_smallest(x0, x1, fx0, fx1)
    quad_ctr = one(x1)
    init_state!(state, b, a, T[quad_ctr], fb, fa, S[])
    state.fnevals = 2
    nothing
end



function init_options(::Order0,
                      state::UnivariateZeroState{T,S};
                      maxevals = nothing,
                      kwargs...) where {T, S}

    # same as order 1 save maxevals
    options = init_options(Order1(), state; kwargs...)
    options.maxevals = maxevals == nothing ? 5 * ceil(Int, -log(eps(T))) : maxevals
    options
end


function assess_convergence(method::Order0, state::UnivariateZeroState{T,S}, options) where {T,S}

    quad_ctr = state.m[1] #add in test on quad qtr
    if quad_ctr > 3
        state.stopped = true
        return true
    end

    # how to call super?
    assess_convergence(Order1(), state, options)

end


function run_bisection(f, xs, state, options)
    #M = Bisection()  # exact for floating point
    M = AlefeldPotraShi() # *usually* exact
    #M = Brent()          # a bit faster, but not always convergent, as implemented (cf. RootTesting)
    steps = state.steps
    init_state!(state, M, f, xs); state.steps += steps
    init_options!(options, M)
    find_zero(M, f, options, state)
    a, b = xs
    u,v = a > b ? (b, a) : (a, b)
    state.message = "Bisection used over ($u, $v), steps not shown"
    return nothing
end

# main algorithm
function update_state(method::Order0, f, state::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T, b::T, quad_ctr = state.xn0, state.xn1, state.m[1]
    fa::S, fb::S = state.fxn0, state.fxn1

    ## we keep fb < fa Could do in init_state to save a conditional XXX
    if abs(fa) < abs(fb)
        a, b= b, a
        fa,fb = fb, fa
    end

    isbracket(fa, fb) && return run_bisection(f, (a, b), state, options)
    
    gamma::T = secant_step(a,b,fa,fb)
    # modify if gamma is too small or too big
        if isinf(gamma) || isnan(gamma) || iszero(abs(gamma-b))
            gamma = b + sign(gamma-b) * 1/1000 * abs(b-a)  # too small
        elseif abs(gamma-b)  >= 100 * abs(b-a)
            gamma = b + sign(gamma-b) * 100 * abs(b-a)  ## too big
        end
    fgamma::T = f(gamma)
    incfn(state)

    
    isbracket(fgamma, fb) && return run_bisection(f, (gamma, b), state, options)


    # decreasing is good
    if abs(fgamma) < abs(fb)
        state.xn0, state.xn1 = b, gamma
        state.fxn0, state.fxn1 = fb, fgamma
        return nothing
    end

    
    # try quad step
    gamma = quad_vertex(a,fa,b, fb, gamma, fgamma)
    fgamma = f(gamma)
    incfn(state)
    state.m[1] += 1 # increment quad_ctr

    if abs(fgamma) < abs(fb)
        state.xn0, state.xn1 = b, gamma
        state.fxn0, state.fxn1 = fb, fgamma
    else
        state.xn0 = gamma
        state.fxn0 = fgamma
    end

    return nothing    
end


##################################################

## Secant
## https://en.wikipedia.org/wiki/Secant_method
""" 
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


