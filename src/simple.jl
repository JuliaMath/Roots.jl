# some simpler (and faster) implementations for root finding
#
# Not exported
#
# These avoid the setup costs of the `find_zero` method, so should be faster
# though they will take similar number of function calls.
#
# `Roots.bisection(f, a, b)`  (Bisection).
# `Roots.secant_method(f, xs)` (Order1) secant method
# `Roots.dfree(f, xs)`  (Order0) more robust secant method
#

## Bisection
##
## Essentially from Jason Merrill https://gist.github.com/jwmerrill/9012954
## cf. http://squishythinking.com/2014/02/22/bisecting-floats/
## This also borrows a trick from https://discourse.julialang.org/t/simple-and-fast-bisection/14886
## where we keep x1 so that y1 is negative, and x2 so that y2 is positive
## this allows the use of signbit over y1*y2 < 0 which avoid < and a multiplication
## this has a small, but noticeable impact on performance.
"""
    bisection(f, a, b; [xatol, xrtol])

Performs bisection method to find a zero of a continuous
function.

It is assumed that `(a,b)` is a bracket, that is, the function has
different signs at `a` and `b`. The interval `(a,b)` is converted to floating point
and shrunk when `a` or `b` is infinite. The function `f` may be infinite for
the typical case. If `f` is not continuous, the algorithm may find
jumping points over the x axis, not just zeros.


If non-trivial tolerances are specified, the process will terminate
when the bracket `(a,b)` satisfies `isapprox(a, b, atol=xatol,
rtol=xrtol)`. For zero tolerances, the default, for `Float64`, `Float32`,
or `Float16` values, the process will terminate at a value `x` with
`f(x)=0` or `f(x)*f(prevfloat(x)) < 0 ` or `f(x) * f(nextfloat(x)) <
0`. For other number types, the `Roots.A42` method is used.

"""
function bisection(f, a::Number, b::Number; xatol=nothing, xrtol=nothing)
    x1, x2 = adjust_bracket(float.((a, b)))
    T = eltype(x1)

    atol = xatol === nothing ? zero(T) : abs(xatol)
    rtol = xrtol === nothing ? zero(one(T)) : abs(xrtol)
    CT = iszero(atol) && iszero(rtol) ? Val(:exact) : Val(:inexact)

    x1, x2 = float(x1), float(x2)
    y1, y2 = f(x1), f(x2)

    _unitless(y1 * y2) >= 0 && error("the interval provided does not bracket a root")

    if isneg(y2)
        x1, x2, y1, y2 = x2, x1, y2, y1
    end

    xm = Roots._middle(x1, x2) # for possibly mixed sign x1, x2
    ym = f(xm)

    while true
        if has_converged(CT, x1, x2, xm, ym, atol, rtol)
            return xm
        end

        if isneg(ym)
            x1, y1 = xm, ym
        else
            x2, y2 = xm, ym
        end

        xm = Roots.__middle(x1, x2)
        ym = f(xm)
    end
end

# -0.0 not returned by __middle, so isneg true on [-Inf, 0.0)
@inline isneg(x::T) where {T<:AbstractFloat} = signbit(x)
@inline isneg(x) = _unitless(x) < 0

@inline function has_converged(::Val{:exact}, x1, x2, m, ym, atol, rtol)
    iszero(ym) && return true
    isnan(ym) && return true
    x1 != m && m != x2 && return false
    return true
end

@inline function has_converged(::Val{:inexact}, x1, x2, m, ym, atol, rtol)
    iszero(ym) && return true
    isnan(ym) && return true
    val = abs(x1 - x2) <= atol + max(abs(x1), abs(x2)) * rtol
    return val
end

"""
    a42(f, ab; atol=nothing, rtol=nothing, λ=0.7, μ = 0.5)

Direct implemenation of Alefeld, Potra, and Shi's Algorithm 4.2. See also [`A42()`](@ref).

* `f`: function to find zero of. (If `f` is 4-times continuously differentiable, convergence to a simple root will be like ``(2 + 7^{1/2})^{1/3} = 1.6686...``
* `ab`: a *bracketing interval
* `atol`, `rtol`: optional tolerances. These are `0` and `eps` respectively by default.


"""
function a42(f, ab; atol=nothing, rtol=nothing, λ=0.7, μ=0.5)
    a, b = adjust_bracket(ab)
    δ₀ = b - a
    fa, fb = f(a), f(b)
    assert_bracket(fa, fb)

    tols = (
        λ    = λ,
        atol = isnothing(atol) ? zero(one(a)) : atol,
        rtol = isnothing(rtol) ? eps(one(a)) : rtol,
    )
    c = a - fa * (b - a) / (fb - fa)
    c = avoid_boundaries(a, c, b, fa, fb, tols)

    fc = f(c)
    iszero(fc) && return c
    e, fee = c, fc
    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

    n = 2
    while true
        δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
        (b - a) ≤ δ && return (abs(fa) < abs(fb) ? a : b)

        ee, fee = d, fd
        for k in 1:2
            if n == 2 || iszero(_pairwise_prod(fa, fb, fd, fee))
                c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
            else
                c = ipzero(a, b, d, ee, fa, fb, fd, fee)
                if (c <= a || b <= c)
                    c = newton_quadratic(a, b, d, fa, fb, fd, k + 1)
                end
            end
            n += 1
            c = avoid_boundaries(a, c, b, fa, fb, tols)
            fc = f(c)
            iszero(fc) && return c

            ee, fee = d, fd
            a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)

            δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
            (b - a) ≤ 2δ && return (abs(fa) < abs(fb) ? a : b)
        end

        n += 1

        u, fu = abs(fa) < abs(fb) ? (a, fa) : (b, fb)
        c = u - 2 * fu * (b - a) / (fb - fa)

        if 2abs(c - u) > (b - a)
            c = a / 2 + b / 2
        end

        c = avoid_boundaries(a, c, b, fa, fb, tols)
        fc = f(c)
        iszero(fc) && return c

        ee, fee = d, fd
        a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
        δ = tolₑ(a, b, fa, fb, tols.atol, tols.rtol)
        (b - a) ≤ 2δ && return (abs(fa) < abs(fb) ? a : b)

        if (b - a) ≥ μ * δ₀
            c = a / 2 + b / 2
            fc = f(c)
            iszero(fc) && return c
            ee, fee = d, fd
            a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb, fc)
        end
        n += 1
    end
end

"""
    secant_method(f, xs; [atol=0.0, rtol=8eps(), maxevals=1000])

Perform secant method to solve `f(x) = 0.`

The secant method is an iterative method with update step
given by `b - fb/m` where `m` is the slope of the secant line between
`(a,fa)` and `(b,fb)`.

The initial values can be specified as a pair of 2, as in `(x₀, x₁)` or
`[x₀, x₁]`, or as a single value, `x₁` in which case a value of `x₀` is chosen.

The algorithm returns m when `abs(fm) <= max(atol, abs(m) * rtol)`.
If this doesn't occur before `maxevals` steps or the algorithm
encounters an issue, a value of `NaN` is returned. If too many steps
are taken, the current value is checked to see if there is a sign
change for neighboring floating point values.

The `Order1` method for `find_zero` also implements the secant
method. This one should be slightly faster, as there are fewer setup costs.

Examples:

```julia
Roots.secant_method(sin, (3,4))
Roots.secant_method(x -> x^5 -x - 1, 1.1)
```

!!! note "Specialization"
    This function will specialize on the function `f`, so that the initial
    call can take more time than a call to the `Order1()` method, though
    subsequent calls will be much faster.  Using `FunctionWrappers.jl` can
    ensure that the initial call is also equally as fast as subsequent
    ones.

"""
function secant_method(
    f,
    xs;
    atol=zero(float(real(first(xs)))),
    rtol=8eps(one(float(real(first(xs))))),
    maxevals=100,
)
    if length(xs) == 1 # secant needs x0, x1; only x0 given
        a = float(xs[1])

        h = eps(one(real(a)))^(1 / 3)
        da = h * oneunit(a) + abs(a) * h^2 # adjust for if eps(a) > h
        b = a + da

    else
        a, b = promote(float(xs[1]), float(xs[2]))
    end
    secant(f, a, b, atol, rtol, maxevals)
end

function secant(f, a::T, b::T, atol=zero(T), rtol=8eps(T), maxevals=100) where {T}
    nan = (0a) / (0a)
    cnt = 0

    fa, fb = f(a), f(b)
    fb == fa && return nan

    uatol = atol / oneunit(atol) * oneunit(real(a))
    adjustunit = oneunit(real(fb)) / oneunit(real(b))

    while cnt < maxevals
        m = b - (b - a) * fb / (fb - fa)
        fm = f(m)

        iszero(fm) && return m
        isnan(fm) && return nan
        abs(fm) <= adjustunit * max(uatol, abs(m) * rtol) && return m
        if fm == fb
            sign(fm) * sign(f(nextfloat(m))) <= 0 && return m
            sign(fm) * sign(f(prevfloat(m))) <= 0 && return m
            return nan
        end

        a, b, fa, fb = b, m, fb, fm

        cnt += 1
    end

    return nan
end

"""
    muller(f, xᵢ; xatol=nothing, xrtol=nothing, maxevals=100)
    muller(f, xᵢ₋₂, xᵢ₋₁, xᵢ; xatol=nothing, xrtol=nothing, maxevals=100)

> *Muller’s method* generalizes the secant method, but uses quadratic
> interpolation among three points instead of linear interpolation between two.
> Solving for the zeros of the quadratic allows the method to find complex
> pairs of roots.
> Given three previous guesses for the root `xᵢ₋₂`, `xᵢ₋₁`, `xᵢ`, and the values
> of the polynomial `f` at those points, the next approximation `xᵢ₊₁` is produced.

Excerpt and the algorithm taken from

> W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery
> *Numerical Recipes in C*, Cambridge University Press (2002), p. 371

Convergence here is decided by `xᵢ₊₁ ≈ xᵢ` using the tolerances specified,
which both default to `eps(one(typeof(abs(xᵢ))))^4/5` in the appropriate units.
Each iteration performs three evaluations of `f`.
The first method picks two remaining points at random in relative proximity of `xᵢ`.

Note that the method may return complex result even for real initial values
as this depends on the function.

Examples:
```
muller(x->x^3-1, 0.5, 0.5im, -0.5) # → -0.500 + 0.866…im
muller(x->x^2+2, 0.0, 0.5, 1.0) # → ≈ 0.00 - 1.41…im
muller(x->(x-5)*x*(x+5), rand(3)...) # → ≈ 0.00
muller(x->x^3-1, 1.5, 1.0, 2.0) # → 2.0, Not converged
```

"""
muller(f, x₀::T; kwargs...) where {T} = muller(f, (rand(T, 2) .* x₀)..., x₀; kwargs...)

muller(f, xᵢ₋₂, xᵢ₋₁, xᵢ; kwargs...) = muller(f, promote(xᵢ₋₂, xᵢ₋₁, xᵢ)...; kwargs...)

function muller(
    f,
    oldest::T,
    older::T,
    old::T;
    xatol=nothing,
    xrtol=nothing,
    maxevals=300,
) where {T}
    @assert old ≠ older ≠ oldest ≠ old # we want q to be non-degenerate
    xᵢ₋₂, xᵢ₋₁, xᵢ = oldest, older, old
    fxᵢ₋₂, fxᵢ₋₁ = f(xᵢ₋₂), f(xᵢ₋₁)

    RT = typeof(abs(oldest))
    atol = xatol !== nothing ? xatol : oneunit(RT) * (eps(one(RT)))^(4 / 5)
    rtol = xrtol !== nothing ? xrtol : eps(one(RT))^(4 / 5)

    for i in 1:(maxevals ÷ 3)
        # one evaluation per iteration
        fxᵢ = f(xᵢ)
        x = muller_step(xᵢ₋₂, xᵢ₋₁, xᵢ, fxᵢ₋₂, fxᵢ₋₁, fxᵢ)

        if isnan(x)
            @warn "The algorithm might not have converged, stopping at i=$i:" abs(xᵢ - xᵢ₋₁)
            return xᵢ
        end

        # @debug "Iteration $i:" xᵢ₋₂ xᵢ₋₁ xᵢ x abs(x-xᵢ)
        xᵢ₋₂, xᵢ₋₁, xᵢ = xᵢ₋₁, xᵢ, x
        fxᵢ₋₂, fxᵢ₋₁ = fxᵢ₋₁, fxᵢ
        #stopping criterion
        isapprox(xᵢ, xᵢ₋₁, atol=atol, rtol=rtol) && return xᵢ
    end
    @warn "The algorithm might not have converged, maxevals limit hit:" abs(xᵢ₋₁ - xᵢ)
    return xᵢ
end

function muller_step(a, b, c, fa, fb, fc)
    a, b, c = promote(a, b, c)
    q = qq(a, b, c)
    q² = q^2
    q1 = q + one(q)

    A = q * fc - q * q1 * fb + q² * fa
    B = (q1 + q) * fc - q1^2 * fb + q² * fa
    C = q1 * fc

    den = let
        Δ = B^2 - 4A * C
        typeof(Δ) <: Real &&
            Δ < 0 &&
            throw(
                DomainError(
                    Δ,
                    "Discriminant is negative and the function most likely has complex roots. You might want to call muller with complex input.",
                ),
            )
        Δ = √Δ
        d⁺ = B + Δ
        d⁻ = B - Δ
        abs(d⁺) > abs(d⁻) ? d⁺ : d⁻
    end
    return c - (c - b) * 2C / den
end

@inline qq(a, b, c) = (c - b) / (b - a)

"""
    newton((f, f'), x0; xatol=nothing, xrtol=nothing, maxevals=100)
    newton(fΔf, x0; xatol=nothing, xrtol=nothing, maxevals=100)

Newton's method.

Function may be passed in as a tuple (f, f') *or* as function which returns (f,f/f').

Examples:
```
newton((sin, cos), 3.0)
newton(x -> (sin(x), sin(x)/cos(x)), 3.0, xatol=1e-10, xrtol=1e-10)
```

Note: unlike the call `newton(f, fp, x0)`--which dispatches to a method of `find_zero`, these
two interfaces will specialize on the function that is passed in. This means, these functions
will be faster for subsequent calls, but may be slower for an initial call.

Convergence here is decided by x_n ≈ x_{n-1} using the tolerances specified, which both default to
`eps(T)^4/5` in the appropriate units.

If the convergence fails, will return a `ConvergenceFailed` error.

"""
struct TupleWrapper{F,Fp}
    f::F
    fp::Fp
end
(F::TupleWrapper)(x) = begin
    u, v = F.f(x), F.fp(x)
    return (u, u / v)
end

newton(f::Tuple, x0; kwargs...) = newton(TupleWrapper(f[1], f[2]), x0; kwargs...)
function newton(f, x0; xatol=nothing, xrtol=nothing, maxevals=100)
    x = float(x0)
    T = typeof(x)
    atol = xatol !== nothing ? xatol : oneunit(T) * (eps(one(T)))^(4 / 5)
    rtol = xrtol !== nothing ? xrtol : eps(one(T))^(4 / 5)

    xo = Inf
    for i in 1:maxevals
        fx, Δx = f(x)
        iszero(fx) && return x

        x -= Δx

        if isapprox(x, xo, atol=atol, rtol=rtol)
            return x
        end

        xo = x
    end

    throw(ConvergenceFailed("No convergence"))
end

## --------------------------------------------------

## This is basically Order0(), but with different, default, tolerances employed
## It takes more function calls, but works harder to find exact zeros
## where exact means either iszero(fx), adjacent floats have sign change, or
## abs(fxn) <= 8 eps(xn)
"""
    dfree(f, xs)

A more robust secant method implementation

Solve for `f(x) = 0` using an algorithm from *Personal Calculator Has Key
to Solve Any Equation f(x) = 0*, the SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).

This is also implemented as the `Order0` method for `find_zero`.

The initial values can be specified as a pair of two values, as in
`(a,b)` or `[a,b]`, or as a single value, in which case a value of `b`
is computed, possibly from `fb`.  The basic idea is to follow the
secant method to convergence unless:

* a bracket is found, in which case `AlefeldPotraShi` is used;

* the secant method is not converging, in which case a few steps of a
  quadratic method are used to see if that improves matters.

Convergence occurs when `f(m) == 0`, there is a sign change between
`m` and an adjacent floating point value, or `f(m) <= 2^3*eps(m)`.

A value of `NaN` is returned if the algorithm takes too many steps
before identifying a zero.

# Examples

```julia
Roots.dfree(x -> x^5 - x - 1, 1.0)
```

"""
function dfree(f, xs)
    if length(xs) == 1
        a = float(xs[1])
        fa = f(a)

        h = eps(one(a))^(1 / 3)
        da = h * oneunit(a) + abs(a) * h^2 # adjust for if eps(a) > h
        b = float(a + da)
        fb = f(b)
    else
        a, b = promote(float(xs[1]), float(xs[2]))
        fa, fb = f(a), f(b)
    end

    nan = (0 * a) / (0 * a) # try to preserve type
    cnt, MAXCNT = 0, 5 * ceil(Int, -log(eps(one(a))))  # must be higher for BigFloat
    MAXQUAD = 3

    if abs(fa) > abs(fb)
        a, fa, b, fb = b, fb, a, fa
    end

    # we keep a, b, fa, fb, gamma, fgamma
    quad_ctr = 0
    while !iszero(fb)
        cnt += 1

        if sign(fa) * sign(fb) < 0
            return solve(ZeroProblem(f, (a, b))) # faster than bisection(f, a, b)
        end

        # take a secant step
        gamma = float(b - (b - a) * fb / (fb - fa))
        # modify if gamma is too small or too big
        if iszero(abs(gamma - b))
            gamma = b + 1 / 1000 * abs(b - a)  # too small
        elseif abs(gamma - b) >= 100 * abs(b - a)
            gamma = b + sign(gamma - b) * 100 * abs(b - a)  ## too big
        end
        fgamma = f(gamma)

        # change sign
        if sign(fgamma) * sign(fb) < 0
            return solve(ZeroProblem(f, (gamma, b))) # faster than bisection(f, gamma, b)
        end

        # decreasing
        if abs(fgamma) < abs(fb)
            a, fa, b, fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        gamma = float(quad_vertex(a, fa, b, fb, gamma, fgamma))
        fgamma = f(gamma)
        # decreasing now?
        if abs(fgamma) < abs(fb)
            a, fa, b, fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        quad_ctr += 1
        if (quad_ctr > MAXQUAD) || (cnt > MAXCNT) || iszero(gamma - b) || isnan(gamma)
            bprev, bnext = prevfloat(b), nextfloat(b)
            fbprev, fbnext = f(bprev), f(bnext)
            sign(fb) * sign(fbprev) < 0 && return b
            sign(fb) * sign(fbnext) < 0 && return b
            for (u, fu) in ((b, fb), (bprev, fbprev), (bnext, fbnext))
                abs(fu) / oneunit(fu) <= 2^3 * eps(u / oneunit(u)) && return u
            end
            return nan # Failed.
        end

        if abs(fgamma) < abs(fb)
            b, fb, a, fa = gamma, fgamma, b, fb
        else
            a, fa = gamma, fgamma
        end
    end
    b
end
