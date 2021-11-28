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
@inline function do_guarded_step(
    M::AbstractSecant,
    o::AbstractUnivariateZeroState{T,S},
) where {T,S}
    x, fx = o.xn1, o.fxn1
    1000 * abs(fx) > max(oneunit(S), abs(x) * oneunit(S) / oneunit(T)) * one(T)
end

# check if we should guard against step for method M; call N if yes, P if not
function update_state_guarded(
    M::AbstractSecant,
    N::AbstractUnivariateZeroMethod,
    P::AbstractUnivariateZeroMethod,
    fs,
    o,
    options,
    l=NullTracks(),
)
    if do_guarded_step(M, o)
        return update_state(N, fs, o, options, l)
    else
        update_state(P, fs, o, options, l)
    end
end

##################################################
initial_fncalls(::AbstractSecant) = 2

##################################################

## Secant
## https://en.wikipedia.org/wiki/Secant_method
"""
    Secant()
    Order1()
    Orderφ()


The `Order1()` method is an alias for `Secant`. It specifies the
[secant method](https://en.wikipedia.org/wiki/Secant_method).
This method keeps two values in its state, `xₙ` and `xₙ₋₁`. The
updated point is the intersection point of ``x`` axis with the secant line
formed from the two points. The secant method uses ``1`` function
evaluation per step and has order `φ≈ (1+sqrt(5))/2`.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₂ = f[xᵢ₊₁,xᵢ,α] / f[xᵢ₊₁,xᵢ] * (xᵢ₊₁-α) * (xᵢ - α)`.

"""
struct Secant <: AbstractSecant end
const Order1 = Secant
const Orderφ = Secant

# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractSecant, F::Callable_Function, x)
    x₀, x₁ = x₀x₁(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# initialize from xs, fxs
function init_state(::AbstractSecant, F, x₀, x₁, fx₀, fx₁)
    UnivariateZeroState(x₁, x₀, fx₁, fx₀)
end

function update_state(
    ::Order1,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)

    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1::T = xn1, xn1 - Δ
    fx0, fx1 = fxn1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false
end

### Steffensen
## https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description

"""
    Steffensen()
    Order2()


The quadratically converging
[Steffensen](https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description)
method is used for the derivative-free `Order2()` algorithm. Unlike
the quadratically converging Newton's method, no derivative is
necessary, though like Newton's method, two function calls per step
are. Steffensen's algorithm is more sensitive than Newton's method to
poor initial guesses when `f(x)` is large, due to how `f'(x)` is
approximated. The `Order2` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ - α`, satisfies
`eᵢ₊₁ = f[xᵢ, xᵢ+fᵢ, α] / f[xᵢ,xᵢ+fᵢ] ⋅ (1 - f[xᵢ,α] ⋅ eᵢ²`
"""
struct Steffensen <: AbstractSecant end
struct Order2 <: AbstractSecant end

function update_state(
    method::Order2,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(method, Secant(), Steffensen(), fs, o, options, l)
end

function update_state(
    method::Steffensen,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1

    sgn = sign((fx1 - fx0) / (x1 - x0))
    x2 = x1 - sgn * fx1 / oneunit(S) * oneunit(T)

    f0 = fx1
    f1::S = F(x2)
    incfn(l, 1)

    delta = -sgn * f0 * f0 / (f1 - f0) * oneunit(T) / oneunit(S)

    if isissue(delta)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1::T = x1, x1 - delta
    fx0, fx1 = fx1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false
end

##################################################

"""
    Order5()
    KumarSinghAkanksha()

Implements an order 5 algorithm from *A New Fifth Order Derivative
Free Newton-Type Method for Solving Nonlinear Equations* by Manoj
Kumar, Akhilesh Kumar Singh, and Akanksha, Appl. Math. Inf. Sci. 9,
No. 3, 1507-1513 (2015), DOI: [10.12785/amis/090346](https://doi.org/10.12785/amis/090346). Four function
calls per step are needed.  The `Order5` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ = K₁ ⋅ K₅ ⋅ M ⋅ eᵢ⁵ + O(eᵢ⁶)`

"""
struct Order5 <: AbstractSecant end
struct KumarSinghAkanksha <: AbstractSecant end

function update_state(method::Order5, fs, o::UnivariateZeroState, options, l=NullTracks())
    update_state_guarded(method, Secant(), KumarSinghAkanksha(), fs, o, options, l)
end

function update_state(
    M::KumarSinghAkanksha,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, o.xn1, o.fxn1)

    fwn::S = F(wn)
    incfn(l)

    fp, issue = _fbracket(o.xn1, wn, o.fxn1, fwn)
    if issue
        log_message(l, "Issue with divided difference f[xn, wn]. ")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = wn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fwn

        return o, true
    end

    yn::T = o.xn1 - o.fxn1 / fp
    fyn::S = F(yn)
    incfn(l)

    zn::T = xn - (fxn + fyn) / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_ratio(yn, o.xn1, wn, fyn, o.fxn1, fwn)
    if issue
        log_message(l, "Issue with f[xn,yn]*f[yn,wn] / f[xn, wn]")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = yn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fyn

        return o, true
    end

    @set! o.xn0 = o.xn1
    @set! o.fxn0 = o.fxn1
    x₁::T = zn - fzn / fp
    @set! o.xn1 = x₁
    @set! o.fxn1 = F(o.xn1)
    incfn(l)

    return o, false

    #    nothing
end

struct Order5Derivative <: AbstractSecant end
fn_argout(::Order5Derivative) = 2
function update_state(
    method::Order5Derivative,
    f,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn, fxn = o.xn1, o.fxn1
    a::T, b::S = f(xn)
    fpxn = a / b
    incfn(l)

    if isissue(fpxn)
        return o, true
    end

    yn::T = xn - fxn / fpxn
    fyn::S, Δyn::T = f(yn)
    fpyn = fyn / Δyn
    incfn(l, 2)

    if isissue(fpyn)
        log_message(l, "Issue computing `fpyn`")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = yn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fyn

        return o, true
    end

    zn::T = xn - (fxn + fyn) / fpxn
    fzn::S, _ = f(zn)
    incfn(l, 2)

    xn1::T = zn - fzn / fpyn
    fxn1, _ = f(xn1)
    incfn(l, 2)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1

    return o
end

##################################################

## cf also: https://doi.org/10.1515/tmj-2017-0049
"""
    Order8()
    Thukral8()

Implements an eighth-order algorithm from *New Eighth-Order
Derivative-Free Methods for Solving Nonlinear Equations* by Rajinder
Thukral, International Journal of Mathematics and Mathematical
Sciences Volume 2012 (2012), Article ID 493456, 12 pages DOI:
[10.1155/2012/493456](https://doi.org/10.1155/2012/493456). Four
function calls per step are required.  The `Order8` method replaces a
Steffensen step with a secant step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K ⋅ eᵢ⁸` in
(2.25) of the paper for an explicit `K`.

"""
struct Order8 <: AbstractSecant end
struct Thukral8 <: AbstractSecant end

function update_state(
    method::Order8,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(method, Secant(), Thukral8(), fs, o, options, l)
end

function update_state(
    M::Thukral8,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, xn, fxn)
    fwn::S = F(wn)
    incfn(l)

    if isissue(fwn)
        log_message(l, "issue with Steffensen step fwn")
        @set! o.xn0 = xn
        @set! o.xn1 = wn
        @set! o.fxn0 = fxn
        @set! o.fxn1 = fwn

        return o, true
    end

    fp, issue = _fbracket(xn, wn, fxn, fwn)

    if issue
        log_message(l, "issue with divided difference f[xn, wn]. ")
        return o, true
    end

    yn::T = xn - fxn / fp
    fyn::S = F(yn)
    incfn(l)

    fp, issue = _fbracket(yn, xn, fyn, fxn)
    if issue #fp
        log_message(l, "issue with divided difference f[xn, yn]. ")
        @set! o.xn0 = xn
        @set! o.xn1 = yn
        @set! o.fxn0 = fxn
        @set! o.fxn1 = fyn

        return o, true
    end

    phi = (1 + fyn / fwn)           # pick one of options
    zn::T = yn - phi * fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    if issue
        log_message(l, "issue with divided difference  f[y,z] - f[x,y] + f[x,z]. ")
        @set! o.xn0 = xn
        @set! o.xn1 = zn
        @set! o.fxn0 = fxn
        @set! o.fxn1 = fzn

        return o, true
    end

    w = 1 / (1 - fzn / fwn)

    xi = (1 - 2fyn * fyn * fyn / (fwn * fwn * fxn))

    xn1::T = zn - w * xi * fzn / fp
    fxn1::S = F(xn1)
    incfn(l)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1

    return o, false
end

##################################################

"""
    Order16()
    Thukral16()

Implements the order 16 algorithm from
*New Sixteenth-Order Derivative-Free Methods for Solving Nonlinear Equations*
by R. Thukral,
American Journal of Computational and Applied Mathematics
p-ISSN: 2165-8935;    e-ISSN: 2165-8943; 2012;  2(3): 112-118
DOI: [10.5923/j.ajcam.20120203.08](https://doi.org/10.5923/j.ajcam.20120203.08).

Five function calls per step are required. Though rapidly converging,
this method generally isn't faster (fewer function calls/steps) over
other methods when using `Float64` values, but may be useful for
solving over `BigFloat`.  The `Order16` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K⋅eᵢ¹⁶` for an explicit `K`
in equation (50) of the paper.

"""
struct Order16 <: AbstractSecant end
struct Thukral16 <: AbstractSecant end

function update_state(
    method::Order16,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(method, Secant(), Thukral16(), fs, o, options, l)
end

function update_state(
    M::Thukral16,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1

    wn::T = steff_step(M, xn, fxn)
    fwn::S = F(wn)
    incfn(l)

    fp, issue = _fbracket(xn, wn, fxn, fwn)

    if issue
        log_message(l, "issue with f[xn,wn]")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = wn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fwn
        return o, true
    end

    yn::T = xn - fxn / fp
    fyn::S = F(yn)
    incfn(l)

    fp, issue = _fbracket_ratio(yn, xn, wn, fyn, fxn, fwn)
    if issue
        log_message(l, "issue with f[xn,yn]*f[yn,wn]/f[xn,wn]")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = yn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fyn

        return o, true
    end

    zn::T = yn - fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    u2, u3, u4 = fzn / fwn, fyn / fxn, fyn / fwn

    eta = 1 / (1 + 2 * u3 * u4^2) / (1 - u2)
    if issue
        log_message(l, "Approximate derivative failed")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = zn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fzn

        return o, true
    end

    an::T = zn - eta * fzn / fp
    fan::S = F(an)
    incfn(l)

    fp, issue = _fbracket_ratio(an, yn, zn, fan, fyn, fzn)
    if issue
        log_message(l, "Approximate derivative failed")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = an
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fan

        return o, true
    end

    u1, u5, u6 = fzn / fxn, fan / fxn, fan / fwn

    fp1, issue = _fbracket(xn, yn, fxn, fyn)

    sigma =
        1 + u1 * u2 - u1 * u3 * u4^2 +
        u5 +
        u6 +
        u1^2 * u4 +
        u2^2 * u3 +
        3 * u1 * u4^2 * (u3^2 - u4^2) / (fp1 / oneunit(fp1))

    xn1::T = an - sigma * fan / fp
    fxn1::S = F(xn1)
    incfn(l)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1

    return o, false
end

##################################################
## some means of guarding against large fx when taking a secant step
## TODO: rework this
function steff_step(M::Union{Order5,Order8,Order16}, x::S, fx::T) where {S,T}
    xbar, fxbar = real(x / oneunit(x)), fx / oneunit(fx)
    thresh = max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6

    out = abs(fxbar) <= thresh ? fxbar : sign(fx) * thresh
    x + out * oneunit(x)
end
