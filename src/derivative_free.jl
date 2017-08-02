# Many derivative free methods of different orders

##################################################
## Helpers for the various methods
## issue with approx derivative
isissue(x) = (x == 0.0) || isnan(x) || isinf(x)


"""
heuristic to get a decent first step with Steffensen steps
"""
function steff_step{T}(x::T, fx)
    thresh =  max(1, norm(x)) * sqrt(eps(T)) # max(1, sqrt(abs(x/fx))) * 1e-6
    norm(fx) <= thresh ? fx : sign(fx) * thresh
end

function guarded_secant_step(alpha, beta, falpha, fbeta)
    fp = (fbeta - falpha) /  (beta - alpha)
    Δ = fbeta / fp

    if norm(Δ) >= 100 * norm(alpha - beta) # guard runaway
        Δ = sign(Δ) * 100 * norm(alpha - beta)
    end

    beta - Δ, isissue(Δ)

end


## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    num, den = fb-fa, b - a
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
type Order0 <: AbstractSecant end
type Secant <: AbstractSecant end
const Order1 = Secant

function init_state{T <: AbstractFloat}(method::AbstractSecant, fs, x::Union{T, Vector{T}}, bracket)

    if isa(x, Vector)
        x0, x1 = x[1:2]
        fx0, fx1 = fs.f(x0), fs.f(x1)
    else
        x0 = float(x)
        fx0 = fs.f(x0)
        stepsize = max(1/100, min(abs(fx0), abs(x0/100)))
        x1 = x0 + stepsize
        x0, x1, fx0, fx1  = x1, x0, fs.f(x1), fx0 # switch
    end

    state = UnivariateZeroStateBase(
                                    promote(float(x1), float(x0))...,
                                    promote(fx1, fx0)...,
                                    isa(bracket, Nullable) ? bracket : Nullable(convert(Vector{T}, sort(bracket))),
                                    0, 2,
                                    false, false, false, false,
                                    "")
    state
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
function update_state{T}(method::Order0, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)

    f = fs.f
    alpha, beta = o.xn0, o.xn1
    falpha, fbeta = o.fxn0, o.fxn1
    S = eltype(falpha)

    incsteps(o)

    if sign(falpha) * sign(fbeta) < 0.0
        # use bisection
        verbose = options.verbose; options.verbose=false # turn off verbose
        find_zero(Bisection(), fs, o, options)
        options.verbose = verbose
        o.message = "Used bisection for last step: [a,b] = [$alpha, $beta]"
        return nothing
    end

    gamma, issue = guarded_secant_step(alpha, beta, falpha, fbeta)
    if issue
        o.message = "error with guarded secant step"
        o.stopped = true
        return nothing
    end

    fgamma = f(gamma); incfn(o)
    if sign(fgamma)*sign(fbeta) < 0.0
        o.xn0, o.xn1 = gamma, beta
        o.fxn0, o.fxn1 = fgamma, fbeta
        verbose = options.verbose; options.verbose=false # turn off verbose
        find_zero(Bisection(), fs, o, options)
        options.verbose=verbose
        o.message = "Used bisection for last step"
        return nothing
    end

    if norm(fgamma) <= norm(fbeta)
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


        fgamma = f(gamma); incfn(o)
        incfn(o)

        if norm(fgamma) < norm(fbeta)
            o.xn0, o.xn1 = beta, gamma
            o.fxn0, o.fxn1 = fbeta, fgamma
            return nothing
        end

        theta, issue = guarded_secant_step(beta, gamma, fbeta, fgamma)
        ftheta = f(theta); incfn(o)

        if sign(ftheta) * sign(fbeta) < 0
            o.xn0, o.xn1 = beta, theta
            o.fxn0, o.fxn1 = fbeta, ftheta

            opts = deepcopy(options); #opts.verbose=false
            o.message = "used bisection for last step"
            find_zero(Bisection(), fs, o, opts)
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
function update_state{T}(method::Secant, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    S = eltype(fxn0)

    incsteps(o)

    fp, issue = _fbracket(xn0, xn1, fxn0, fxn1)
    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end

    xn2::T = xn1 -  fxn1 / fp
    fxn2::S = fs.f(xn2)
    incfn(o)

    o.xn0, o.xn1 = xn1, xn2
    o.fxn0, o.fxn1 = fxn1, fxn2


    nothing
end

"""
secant_method: solve for zero of `f(x) = 0` using the secant method
"""
secant_method(f, x0::Real, x1::Real; kwargs...) = find_zero(f, map(float, [x0,x1]), Order1(); kwargs...)


##################################################

### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description
type Steffensen <: UnivariateZeroMethod
end
const Order2 = Steffensen

function update_state{T}(method::Steffensen, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)

    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
    incfn(o)

    fp, issue = _fbracket(xn, wn, fxn, fwn)

    if issue
        o.stopped = true
        o.message = "Derivative approximation had issues"
        return
    end

    xn1::T = xn - fxn / fp


    fxn1::S = fs.f(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1


    nothing
end

steffenson(f, x0; kwargs...) = find_zero(f, x0, Steffensen(); kwargs...)

##################################################

## A New Fifth Order Derivative Free Newton-Type Method for Solving Nonlinear Equations
## Manoj Kumar, Akhilesh Kumar Singh, and Akanksha Srivastava
## Appl. Math. Inf. Sci. 9, No. 3, 1507-1513 (2015)
## http://www.naturalspublishing.com/files/published/ahb21733nf19a5.pdf
type Order5 <: UnivariateZeroMethod
end

function update_state{T}(method::Order5, fs::DerivativeFree, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
    incfn(o)

    fp, issue = _fbracket(xn, wn, fxn, fwn)
    if issue
        o.xn0, o.xn1 = xn, wn
        o.fxn0, o.fxn1 = fxn, fwn
        o.message = "Issue with divided difference f[xn, wn]"
        o.stopped  = true
        return
    end

    yn::T = xn - fxn / fp
    fyn::S = fs.f(yn)
    incfn(o)


    zn::T = xn - (fxn + fyn) / fp
    fzn::S = fs.f(zn)
    incfn(o)

    fp, issue = _fbracket_ratio(yn, xn, wn, fyn, fxn, fwn)
    if issue
        o.xn0, o.xn1 = xn, yn
        o.fxn0, o.fxn1 = fxn, fyn
        o.message = "Issue with f[xn,yn]*f[yn,wn] / f[xn, wn]"
        o.stopped = true
        return
    end

    xn1::T = zn  - fzn  / fp
    fxn1::S = fs.f(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end

## If we have a derivative
function update_state{T}(method::Order5, fs::FirstDerivative, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)


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


# kss5(f, x0; kwargs...) = kss5(f, D(f), x0; kwargs...)
# function kss5(f, fp, x0; kwargs...)
#     fs = FirstDerivative(f, fp)
#     T = typeof(float(x0))
#     state = init_state(Order5(), fs, x0, Nullable{Vector{T}}())
#     options = UnivariateZeroOptions(eps(T), eps(T), eps(T), eps(T), 40, 200, true)
#     zp = UnivariateZeroProblem(fs, state, options)

#     find_zero(zp, Order5())
# end


##################################################

## New Eighth-Order Derivative-Free Methods for Solving Nonlinear Equations
## Rajinder Thukral
## International Journal of Mathematics and Mathematical Sciences
## Volume 2012 (2012), Article ID 493456, 12 pages
## http://dx.doi.org/10.1155/2012/493456

type Order8 <: UnivariateZeroMethod
end

function update_state{T}(method::Order8, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
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
    fyn::S = fs.f(yn)
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
    fzn::S = fs.f(zn)
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
    fxn1::S = fs.f(xn1)
    incfn(o)

    o.xn0,o.xn1 = xn, xn1
    o.fxn0,o.fxn1 = fxn, fxn1

    nothing
end

##################################################

## New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations
## R. Thukral
## American Journal of Computational and Applied Mathematics
## p-ISSN: 2165-8935    e-ISSN: 2165-8943
## 2012;  2(3): 112-118
## doi: 10.5923/j.ajcam.20120203.08

type Order16 <: UnivariateZeroMethod
end

function update_state{T}(method::Order16, fs, o::UnivariateZeroState{T}, options::UnivariateZeroOptions)
    xn = o.xn1
    fxn = o.fxn1
    S = eltype(fxn)

    incsteps(o)

    wn::T = xn + steff_step(xn, fxn)
    fwn::S = fs.f(wn)
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
    fyn::S = fs.f(yn)
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
    fzn::S = fs.f(zn)
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
    fan = fs.f(an)
    incfn(o)

    fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn
    sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
        u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/_fbracket(xn,yn, fxn, fyn)[1]

    if issue
        o.xn0, o.xn1 = xn, an
        o.fxn0, o.fxn1 = fxn, fan
        o.stopped = true
        o.message = "Approximate derivative failed"
        return
    end

    xn1 = an - sigma * fan / fp
    fxn1 = fs.f(xn1)
    incfn(o)

    o.xn0, o.xn1 = xn, xn1
    o.fxn0, o.fxn1 = fxn, fxn1

    nothing
end


