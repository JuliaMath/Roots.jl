# Many derivative free methods of different orders
#
# TODO: rework Order5 #https://pdfs.semanticscholar.org/ce50/3210d96f653a14b28da96600d5990d2abe97.pdf
# https://content.sciendo.com/view/journals/tmj/10/4/article-p103.xml 7 and 8
# order8: https://www.hindawi.com/journals/ijmms/2012/493456/ref/



##################################################
## Guard against non-robust algorithms
##
## By default, we do this by deciding if we
## should take a secant step instead of the algorithm For example, for
## Steffensen which is quadratically convergent and the Secant method
## which is only superlinear,
## the error, e_{n+1} = x_{n+1} - alpha, may be smaller after a secant
## step than a Steffensen step. (It is only once x_n is close enough
## to alpha that the method is quadratically convergent.
## The Steffensen error is
## Δn+1 = f[x,x+fx, alpha]/f[x, x+fx] * (1 - f[x, alpha]) (x-alpha)^2
##      ≈ f''/(2f') * ( 1 + f') Δn^2
## The Secant error is
## Δn+1 = f[x,x_{-1},alpha] / f[x,x_{-1}] * (x-alpha) * (x_{-1} - alpha)
##      ≈  f''/(2f')  Δn ⋅ Δn-1
## The ratio is ≈ (1 + f')(Δn / Δn-1)
## It seems reasonable, that a Steffensen step is preferred when
## the ratio satisfies -1 < (1+f') ⋅ Δn /Δn-1 < 1
## We could use f' ~ fp = (fx1-fx0)/(x1-x0); but our proxy for
## Δn/Δn-1 is problematic, as we don't know alpha, and using xn-x_{n-1}
## can be an issue when only x1 and not x0 is specified. This needs
## working around.
##
## Instead, as Steffensen is related to Newton as much as
## (f(x+fx) - fx)/fx  ≈ f'(x), we take a Steffensen step if |fx|
## is small enough. For this we use |fx| <= x/1000; which
## seems to work reasonably well over several different test cases.
@inline function do_guarded_step(M::AbstractSecant, o::UnivariateZeroState{T,S}) where {T, S}
    x, fx = o.xn1, o.fxn1
    1000 * abs(fx) >  max(oneunit(S), abs(x) * oneunit(S) /oneunit(T)) * one(T)
end


# check if we should guard against step for method M; call N if yes, P if not
function update_state_guarded(M::AbstractSecant,N::AbstractUnivariateZeroMethod, P::AbstractUnivariateZeroMethod, fs, o, options)
    if do_guarded_step(M, o)
        #@debug "do secant step"
        return update_state(N, fs, o, options)
    else
        #@debug "do $N step"
        update_state(P, fs, o, options)
    end
end

##################################################

function init_state(method::AbstractSecant, fs, x::Number)
    x1 = float(x)
    x0 = _default_secant_step(x1)
    init_state(method, fs, (x0, x1))
end

function init_state(method::AbstractSecant, fs, x::Union{Tuple, Vector})
    x0, x1 = promote(float(x[1]), float(x[2]))
    fx0, fx1 = fs(x0), fs(x1)
    state = UnivariateZeroState(x1, x0, zero(x1)/zero(x1)*oneunit(x1), eltype(x1)[],
                                fx1, fx0, fx1, eltype(fx1)[],
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
slower, alternative to the other derivative-free root-finding
methods. The implementation roughly follows the algorithm described in
*Personal Calculator Has Key to Solve Any Equation f(x) = 0*, the
SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
The basic idea is to use a secant step. If along the way a bracket is
found, switch to bisection, using `AlefeldPotraShi`.  If the secant
step fails to decrease the function value, a quadratic step is used up
to 4 times.

This is not really 0-order: the secant method has order
1.6...[https://en.wikipedia.org/wiki/Secant_method#Comparison_with_other_root-finding_methods]
and the the bracketing method has order
1.6180...[http://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192965-2/S0025-5718-1993-1192965-2.pdf]
so for reasonable starting points, this algorithm should be
superlinear, and relatively robust to non-reasonable starting points.

"""
struct Order0 <: AbstractSecant end

function find_zero(fs, x0, method::Order0;
                   tracks::AbstractTracks=NullTracks(),
                   verbose=false,
                   kwargs...)
    M = Order1()
    N = AlefeldPotraShi()

    find_zero(callable_function(fs), x0, M, N; tracks=tracks,verbose=verbose, kwargs...)
end

##################################################

## Secant
## https://en.wikipedia.org/wiki/Secant_method
"""
    Order1()
    Roots.Secant()

The `Order1()` method is an alias for `Secant`. It specifies the
[secant method](https://en.wikipedia.org/wiki/Secant_method).
This method keeps two values in its state, `x_n` and `x_n1`. The
updated point is the intersection point of x axis with the secant line
formed from the two points. The secant method uses 1 function
evaluation per step and has order `(1+sqrt(5))/2`.

The error, `e_n = x_n - alpha`, satisfies
`e2 = f[x1,x0,alpha] / f[x1,x0] * (x1-alpha) * (x0 - alpha)`.

"""
struct Secant <: AbstractSecant end
const Order1 = Secant

function update_state(method::Order1, fs, o::UnivariateZeroState{T,S}, options) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1

    delta = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)


    if isinf(delta) || isnan(delta)
        o.stopped = true
        o.message = "Increment `Δx` has issues. "
         return
     end

    o.xn0 = xn1
    o.xn1 -= delta
    o.fxn0 = fxn1
    o.fxn1 = (tmp::S = fs(o.xn1))
    incfn(o)

    nothing

end

##################################################


"""
    Roots.Order1B()
    Roots.King()

A superlinear (order 1.6...) modification of the secant method for multiple roots.
Presented in A SECANT METHOD FOR MULTIPLE ROOTS, by RICHARD F. KING, BIT 17 (1977), 321-328

The basic idea is similar to Schroder's method: apply the secant method
to  `f/f'`. However, this uses `f' ~ fp = (fx - f(x-fx))/fx` (a Steffensen step). In
this implementation, `Order1B`, when `fx` is too big, a single secant step of `f`
is used.

The *asymptotic* error, `e_n = x_n - alpha`, is given by
`e2 = 1/2⋅G''/G'⋅ e0⋅e1 + (1/6⋅G'''/G' - (1/2⋅G''/G'))^2⋅e0⋅e1⋅(e0+e1)`.

"""
struct Order1B <: AbstractSecant end
struct King <: AbstractSecant end


function update_state(method::Order1B, fs, o::UnivariateZeroState{T,S}, options)  where {T, S}
    update_state_guarded(method, Secant(), King(), fs, o, options)
end


function update_state(method::King, fs,
                      o::UnivariateZeroState{T,S}, options)  where {T, S}


    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1

    # G(x,f0,f_1) = -fx^2/(f_1 - f0)
    f0 = fx1
    f_1::S = fs(x1 - f0*oneunit(T)/oneunit(S))
    incfn(o, 1)

    G1 = -f0^2 / (f_1 - f0)

    if isempty(o.fm)
        f0 = fx0
        tmp::S = fs(x0 - f0*oneunit(T)/oneunit(S))
        f_1 = tmp
        incfn(o,1)
        G0 = -f0^2 / (f_1 - f0)
    else
        G0 = first(o.fm)
    end

    m = (x1-x0)/(G1-G0) # approximate value of `m`, the multiplicity
    if abs(m) <= 1e-2 * oneunit(m)
        #@info "small m estimate, stopping"
        o.stopped  = true
        o.message = "Estimate for multiplicity has issues. "
        return
    end

    delta = G1 * (x1 - x0) / (G1 - G0)
    empty!(o.fm); push!(o.fm, G1)

    if isissue(delta)
        o.stopped  = true
        o.message = "Increment `Δx` has issues. "
        return
    end


    o.xn0, o.fxn0 = o.xn1, o.fxn1
    o.xn1 -= delta
    o.fxn1 = fs(o.xn1)
    incfn(o)

    nothing
end





### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description

"""
    Order2()
    Roots.Steffensen()

The quadratically converging
[Steffensen](https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description)
method is used for the derivative-free `Order2()` algorithm. Unlike
the quadratically converging Newton's method, no derivative is
necessary, though like Newton's method, two function calls per step
are. Steffensen's algorithm is more sensitive than Newton's method to
poor initial guesses when `f(x)` is large, due to how `f'(x)` is
approximated. This algorithm, `Order2`, replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `e_n - alpha`, satisfies
`e1 = f[x0, x+f0, alpha] / f[x0,x0+f0] ⋅ (1 - f[x0,alpha] ⋅ e0^2`
"""
struct Order2 <: AbstractSecant end
struct Steffensen <: AbstractSecant end

function update_state(method::Order2, fs, o::UnivariateZeroState{T,S}, options)  where {T, S}
     update_state_guarded(method, Secant(), Steffensen(), fs, o, options)
end


function update_state(method::Steffensen, fs, o::UnivariateZeroState{T,S}, options)  where {T, S}

    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1

    sgn = sign((fx1 - fx0) / (x1 - x0))
    x2 = x1 - sgn * fx1 / oneunit(S) * oneunit(T)

    f0 = fx1
    f1::S = fs(x2)
    incfn(o, 1)

    delta = -sgn * f0 * f0 / (f1 - f0) * oneunit(T) / oneunit(S)

    if isissue(delta)
        o.stopped  = true
        o.message = "Increment `Δx` has issues. "
        return
    end

    o.xn0, o.fxn0 = o.xn1, o.fxn1
    o.xn1 -= delta
    o.fxn1 = fs(o.xn1)
    incfn(o)

    nothing
end

##################################################

### Order2B() Esser method
"""
    Roots.Order2B()
    Roots.Esser()

Esser's method. This is a quadratically convergent method that, like
Schroder's method, does not depend on the multiplicity of the
zero. Schroder's method has update step `x - r2/(r2-r1) * r1`, where `ri =
f^(i-1)/f^(i)`. Esser approximates `f' ~ f[x-h, x+h], f'' ~
f[x-h,x,x+h]`, where `h = fx`, as with Steffensen's method, Requiring 3
function calls per step. The implementation `Order2B` uses a secant
step when |fx| is considered too large.


Esser, H. Computing (1975) 14: 367. https://doi.org/10.1007/BF02253547
Eine stets quadratisch konvergente Modifikation des Steffensen-Verfahrens


Example
```
f(x) = cos(x) - x
g(x) = f(x)^2
x0 = pi/4
find_zero(f, x0, Order2(), verbose=true)        #  3 steps / 7 function calls
find_zero(f, x0, Roots.Order2B(), verbose=true) #  4 / 9
find_zero(g, x0, Order2(), verbose=true)        #  22 / 45
find_zero(g, x0, Roots.Order2B(), verbose=true) #  4 / 10
```
"""
struct Order2B <: AbstractSecant end
struct Esser <: AbstractSecant end

function update_state(method::Order2B, fs, o::UnivariateZeroState{T,S}, options)  where {T, S}
     update_state_guarded(method, Secant(), Esser(), fs, o, options)
end

function update_state(method::Esser, fs,
                      o::UnivariateZeroState{T,S}, options)  where {T, S}


    x1, fx1 = o.xn1, o.fxn1

    f0 = fx1

    f1  = fs(x1 + f0 * oneunit(T) / oneunit(S))
    f_1 = fs(x1 - f0 * oneunit(T) / oneunit(S))
    incfn(o, 2)

    # h = f0
    # r1 = f/f' ~ f/f[x+h,x-h]
    # r2 = f'/f'' ~ f[x+h, x-h]/f[x-h,x,x+h]
    r1 = f0 * 2*f0 / (f1 - f_1) * oneunit(T) / oneunit(S)
    r2 = (f1 - f_1)/(f1 - 2*f0 + f_1) * f0/2 * oneunit(T) / oneunit(S)

    k = r2/(r2-r1)  # ~ m

    if abs(k) <= 1e-2 * oneunit(k)
        #@info "estimate for m is *too* small"
        o.stopped  = true
        o.message = "Estimate for multiplicity had issues. "
        return
    end

    delta = k * r1

    if isissue(delta)
        o.stopped  = true
        o.message = "Increment `Δx` has issues. "
        return
    end

    o.xn0, o.fxn0 = o.xn1, o.fxn1
    o.xn1 -= delta

    o.fxn1 = fs(o.xn1)
    incfn(o)

    nothing
end




##################################################

"""
    Order5()

Implements an order 5 algorithm from *A New Fifth Order Derivative
Free Newton-Type Method for Solving Nonlinear Equations* by Manoj
Kumar, Akhilesh Kumar Singh, and Akanksha, Appl. Math. Inf. Sci. 9,
No. 3, 1507-1513 (2015), DOI: 10.12785/amis/090346. Four function
calls per step are needed.

The error, `e_n = x_n - alpha`, satisfies
`e1 = K_1 ⋅ K_5 ⋅ M ⋅ e0^5 + O(e0^6)`

"""
struct Order5 <: AbstractSecant end
struct KumarSinghAkanksha <: AbstractSecant end




function update_state(M::Union{Order5, KumarSinghAkanksha}, fs, o::UnivariateZeroState{T,S},
                      options) where {T, S}

    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, o.xn1, o.fxn1)

    fwn::S = fs(wn)
    incfn(o)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)
    if issue
        o.xn0, o.xn1 = o.xn1, wn
        o.fxn0, o.fxn1 = o.fxn1, fwn
        o.message = "Issue with divided difference f[xn, wn]. "
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



## If we have a derivative, we have this. (Deprecate?)
function update_state(method::Order5, fs::Union{FirstDerivative,SecondDerivative},
                      o::UnivariateZeroState{T,S}, options)  where {T, S}


    xn, fxn = o.xn1, o.fxn1
    a::T, b::S = fΔx(fs, xn)
    fpxn = a/b
    incfn(o)

    if isissue(fpxn)
        o.stopped  = true
        return
    end

    yn::T = xn - fxn / fpxn
    fyn::S, Δyn::T = fΔx(fs, yn)
    fpyn = fyn / Δyn
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

##################################################


## cf also: https://doi.org/10.1515/tmj-2017-0049
"""
    Order8()

Implements an eighth-order algorithm from *New Eighth-Order
Derivative-Free Methods for Solving Nonlinear Equations* by Rajinder
Thukral, International Journal of Mathematics and Mathematical
Sciences Volume 2012 (2012), Article ID 493456, 12 pages DOI:
10.1155/2012/493456. Four function calls per step are required.

The error, `e_n = x_n - alpha`, is expressed as `e1 = K ⋅ e0^8` in
(2.25) of the paper for an explicit K.

"""
struct Order8 <: AbstractSecant end
struct Thukral8 <: AbstractSecant end

function update_state(M::Union{Thukral8, Order8}, fs, o::UnivariateZeroState{T,S}, options) where {T, S}

    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, xn, fxn)
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
        o.message = "issue with divided difference f[xn, wn]. "
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
        o.message = "issue with divided difference f[xn, yn]. "
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
        o.message = "issue with divided difference  f[y,z] - f[x,y] + f[x,z]. "
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

Implements the order 16 algorithm from
*New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by R. Thukral,
American Journal of Computational and Applied Mathematics
p-ISSN: 2165-8935;    e-ISSN: 2165-8943; 2012;  2(3): 112-118
doi: 10.5923/j.ajcam.20120203.08.

Five function calls per step are required. Though rapidly converging,
this method generally isn't faster (fewer function calls/steps) over
other methods when using `Float64` values, but may be useful for
solving over `BigFloat`.

The error, `e_n = x_n - alpha`, is expressed as `e1 = K
e_0^16` for an explicit `K` in equation (50) of the paper.

"""
struct Order16 <: AbstractSecant end
struct Thukral16 <: AbstractSecant end

function update_state(M::Union{Thukral16, Order16}, fs, o::UnivariateZeroState{T,S}, options) where {T, S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, xn, fxn)
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

##################################################
## some means of guarding against large fx when taking a secant step
## TODO: rework this
function steff_step(M::Union{Order5, Order8, Order16}, x::S, fx::T) where {S, T}

    xbar, fxbar = real(x/oneunit(x)), fx/oneunit(fx)
    thresh =  max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6

    out = abs(fxbar) <= thresh ? fxbar  : sign(fx) * thresh
    x + out * oneunit(x)

end
