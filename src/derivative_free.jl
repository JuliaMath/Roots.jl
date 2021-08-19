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
@inline function do_guarded_step(M::AbstractSecant, o::AbstractUnivariateZeroState{T,S}) where {T, S}
    x, fx = o.xn1, o.fxn1
    1000 * abs(fx) >  max(oneunit(S), abs(x) * oneunit(S) /oneunit(T)) * one(T)
end


# check if we should guard against step for method M; call N if yes, P if not
function update_state_guarded(M::AbstractSecant,N::AbstractUnivariateZeroMethod, P::AbstractUnivariateZeroMethod, fs, o, options, l=NullTracks())
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
    Order1()
    Roots.Secant()

The `Order1()` method is an alias for `Secant`. It specifies the
[secant method](https://en.wikipedia.org/wiki/Secant_method).
This method keeps two values in its state, `x_n` and `x_n1`. The
updated point is the intersection point of x axis with the secant line
formed from the two points. The secant method uses 1 function
evaluation per step and has order `(1+sqrt(5))/2`.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₂ = f[xᵢ₊₁,xᵢ,α] / f[xᵢ₊₁,xᵢ] * (xᵢ₊₁-α) * (xᵢ - α)`.

"""
struct Secant <: AbstractSecant end
const Order1 = Secant

# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractSecant, F::Callable_Function, x)
    x₀,x₁ = x₀x₁(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# initialize from xs, fxs
function init_state(::AbstractSecant, F, x₀, x₁, fx₀, fx₁)
    UnivariateZeroState(x₁, x₀, fx₁, fx₀)
end



function update_state(method::Order1, F, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks()) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)


    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
     end

    x0,x1 = xn1, xn1 - Δ
    fx0, fx1 = fxn1, F(x1)
    incfn(l)


    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false

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

The *asymptotic* error, `eᵢ = xᵢ - α`, is given by
`eᵢ₊₂ = 1/2⋅G''/G'⋅ eᵢ⋅eᵢ₊₁ + (1/6⋅G'''/G' - (1/2⋅G''/G'))^2⋅eᵢ⋅eᵢ₊₁⋅(eᵢ+eᵢ₊₁)`.

"""
struct King <: AbstractSecant end
struct Order1B <: AbstractSecant end

struct KingState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    G0::S
end


function init_state(M::Union{King, Order1B}, F, x₀, x₁::T, fx₀, fx₁::S) where {T,S}
    fₛ₀ = F(x₀ - fx₀*oneunit(T)/oneunit(S))
    G₀ = -fx₀^2 / (fₛ₀ - fx₀)
    KingState(x₁, x₀, fx₁, fx₀, G₀)
end
initial_fncalls(::Union{King, Order1B}) = 3

function update_state(method::Order1B, F, o::KingState{T,S}, options, l=NullTracks())  where {T, S}
    if do_guarded_step(Order1B(), o)
        state, flag = update_state(Order1(), F, o, options, l)

        x0,fx0 = state.xn0, state.fxn0 # clunky! Need to update G₀ after Order1() step
        fₛ = F(x0 - fx0*oneunit(T)/oneunit(S))
        incfn(l)
        G₀ = -fx0^2 / (fₛ - fx0)
        @set! state.G0 = G₀

        return (state, flag)
    else
        update_state(King(), F, o, options, l)
    end
end


function update_state(method::King, F,
                      o::KingState{T,S}, options, l=NullTracks())  where {T, S}


    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1
    G₀ = o.G0
    fₛ₁ = F(x1 - fx1*oneunit(T)/oneunit(S))
    incfn(l)
    G₁ = -fx1^2 / (fₛ₁ - fx1)

    m = (x1-x0) / (G₁-G₀) # approximate value of `m`, the multiplicity

    if abs(m) <= 1e-2 * oneunit(m)
        log_message(l, "Estimate for multiplicity has issues. ")
        return (o, true)
    end

    Δ = G₁ * (x1 - x0) / (G₁ - G₀)

    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1 = x1, x1 - Δ
    fx0, fx1 = fx1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1
    @set! o.G0 = G₁


    return o, false

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

The error, `eᵢ - α`, satisfies
`eᵢ₊₁ = f[xᵢ, xᵢ+fᵢ, α] / f[xᵢ,xᵢ+fᵢ] ⋅ (1 - f[xᵢ,α] ⋅ eᵢ²`
"""
struct Order2 <: AbstractSecant end
struct Steffensen <: AbstractSecant end

function update_state(method::Order2, fs, o::UnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}
     update_state_guarded(method, Secant(), Steffensen(), fs, o, options, l)
end


function update_state(method::Steffensen, F, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}

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

    x0, x1 = x1, x1 - delta
    fx0, fx1 = fx1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false
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
step when `|fx|` is considered too large.


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
struct Esser <: AbstractSecant end
struct Order2B <: AbstractSecant end

function update_state(method::Order2B, fs, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}
     update_state_guarded(method, Secant(), Esser(), fs, o, options, l)
end

function update_state(method::Esser, F,
                      o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}


    x1, fx1 = o.xn1, o.fxn1

    f0 = fx1

    f1  = F(x1 + f0 * oneunit(T) / oneunit(S))
    f_1 = F(x1 - f0 * oneunit(T) / oneunit(S))
    incfn(l, 2)

    # h = f0
    # r1 = f/f' ~ f/f[x+h,x-h]
    # r2 = f'/f'' ~ f[x+h, x-h]/f[x-h,x,x+h]
    r1 = f0 * 2*f0 / (f1 - f_1) * oneunit(T) / oneunit(S)
    r2 = (f1 - f_1)/(f1 - 2*f0 + f_1) * f0/2 * oneunit(T) / oneunit(S)

    k = r2/(r2-r1)  # ~ m

    if abs(k) <= 1e-2 * oneunit(k)
        log_message(l, "Estimate for multiplicity had issues. ")
        return o, true
#        #@info "estimate for m is *too* small"
    end

    delta = k * r1

    if isissue(delta)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1 = x1, x1 - delta
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

Implements an order 5 algorithm from *A New Fifth Order Derivative
Free Newton-Type Method for Solving Nonlinear Equations* by Manoj
Kumar, Akhilesh Kumar Singh, and Akanksha, Appl. Math. Inf. Sci. 9,
No. 3, 1507-1513 (2015), DOI: [10.12785/amis/090346](https://doi.org/10.12785/amis/090346). Four function
calls per step are needed.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ = K₁ ⋅ K₅ ⋅ M ⋅ eᵢ⁵ + O(eᵢ⁶)`

"""
struct Order5 <: AbstractSecant end
struct KumarSinghAkanksha <: AbstractSecant end




function update_state(M::Union{Order5, KumarSinghAkanksha}, F, o::AbstractUnivariateZeroState{T,S},
                      options, l=NullTracks()) where {T, S}

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
    @set! o.xn1 = zn  - fzn  / fp
    @set! o.fxn1 = F(o.xn1)
    incfn(l)

    return o, false


#    nothing
end

struct Order5Derivative <: AbstractSecant end
fn_argout(::Order5Derivative) = 2
function update_state(method::Order5Derivative, f,
                      o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}


    xn, fxn = o.xn1, o.fxn1
    a::T, b::S = f(xn)
    fpxn = a/b
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


    zn::T = xn  - (fxn + fyn) / fpxn
    fzn::S, _ = f(zn)
    incfn(l, 2)

    xn1 = zn - fzn / fpyn
    fxn1,_ = f(xn1)
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

Implements an eighth-order algorithm from *New Eighth-Order
Derivative-Free Methods for Solving Nonlinear Equations* by Rajinder
Thukral, International Journal of Mathematics and Mathematical
Sciences Volume 2012 (2012), Article ID 493456, 12 pages DOI:
[10.1155/2012/493456](https://doi.org/10.1155/2012/493456). Four function calls per step are required.

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K ⋅ eᵢ⁸` in
(2.25) of the paper for an explicit `K`.

"""
struct Order8 <: AbstractSecant end
struct Thukral8 <: AbstractSecant end

function update_state(M::Union{Thukral8, Order8}, F, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks()) where {T, S}

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
    zn =  yn - phi * fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue =  _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    if issue
        log_message(l, "issue with divided difference  f[y,z] - f[x,y] + f[x,z]. ")
        @set! o.xn0 = xn
        @set! o.xn1 = zn
        @set! o.fxn0 = fxn
        @set! o.fxn1 = fzn

        return o, true
    end

    w = 1 / (1 - fzn/fwn)

    xi = (1 - 2fyn*fyn*fyn / (fwn * fwn * fxn))

    xn1 = zn - w * xi * fzn / fp
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

The error, `eᵢ = xᵢ - α`, is expressed as `eᵢ₊₁ = K⋅eᵢ¹⁶` for an explicit `K`
in equation (50) of the paper.

"""
struct Order16 <: AbstractSecant end
struct Thukral16 <: AbstractSecant end

function update_state(M::Union{Thukral16, Order16}, F, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks()) where {T, S}
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



    zn = yn - fyn / fp
    fzn::S = F(zn)
    incfn(l)

    fp, issue = _fbracket_diff(xn, yn, zn, fxn, fyn, fzn)
    u2, u3, u4 = fzn/fwn, fyn/fxn, fyn/fwn


    eta = 1 / (1 + 2*u3*u4^2) / (1 - u2)
    if issue
        log_message(l, "Approximate derivative failed")
        @set! o.xn0 = o.xn1
        @set! o.xn1 = zn
        @set! o.fxn0 = o.fxn1
        @set! o.fxn1 = fzn

        return o, true
    end

    an = zn - eta * fzn / fp
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

    u1, u5, u6 = fzn/fxn, fan/fxn, fan/fwn

    fp1, issue = _fbracket(xn,yn, fxn, fyn)

    sigma =  1 + u1*u2 - u1*u3*u4^2 + u5 + u6 + u1^2*u4 +
    u2^2*u3 + 3*u1*u4^2*(u3^2 - u4^2)/(fp1/oneunit(fp1))


    xn1 = an - sigma * fan / fp
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
function steff_step(M::Union{Order5, Order8, Order16}, x::S, fx::T) where {S, T}

    xbar, fxbar = real(x/oneunit(x)), fx/oneunit(fx)
    thresh =  max(1, abs(xbar)) * sqrt(eps(one(xbar))) #^(1/2) # max(1, sqrt(abs(x/fx))) * 1e-6

    out = abs(fxbar) <= thresh ? fxbar  : sign(fx) * thresh
    x + out * oneunit(x)

end
