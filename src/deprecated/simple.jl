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


#=
"""
    Roots.a42(f, ab; atol=nothing, rtol=nothing, λ=0.7, μ = 0.5)

Direct implementation of Alefeld, Potra, and Shi's Algorithm 4.2. See also [`A42()`](@ref).

* `f`: function to find zero of. (If `f` is 4-times continuously differentiable, convergence to a simple root will be like ``(2 + 7^{1/2})^{1/3} = 1.6686...``
* `ab`: a *bracketing interval
* `atol`, `rtol`: optional tolerances. These are `0` and `eps` respectively by default.

Not exported
"""
=#
function a42(f, ab; atol=nothing, rtol=nothing, λ=0.7, μ=0.5)
    #Base.depwarn("`a42(f, ab)` is deprecated; use `find_zero(f, ab, Roots.AlefeldPotraShi())` instead.", :a42)
    a, b = adjust_bracket(ab)
    δ₀ = b - a
    fa, fb = f(a), f(b)
    assert_bracket(fa, fb)

    tols = (
        λ   = λ,
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
    modab(f, left, right, args...; kwargs...)

Implementation of "Modified Anderson-Bjork’s method for solving non-linear equations in structural mechanics" by N Ganchovski and A Traykov. Code contributed by @Proektsoft-EOOD in issue [#487](https://github.com/JuliaMath/Roots.jl/issues/487)

This is a modified Anderson-Bjork bracketing algorithm which takes the least number of iterations over a wide-ranging test suite.
"""
function modab(f, left::Real, right::Real, target::Real=0.0; precision::Float64=1e-14, maxIter::Int=200)
    #Base.depwarn("`modab(f, left, right)` is deprecated; use `find_zero(f, (left, right), Roots.ModAB())` instead.", :a42)
    x1, x2 = min(left, right), max(left, right)
    y1 = f(x1) - target
    abs(y1) <= precision && return x1
    y2 = f(x2) - target
    abs(y2) <= precision && return x2
    eps1 = precision * 1e-3
    eps2 = precision * (x2 - x1)
    if abs(target) >= 1
        eps1 *= abs(target)
    else
        eps1 = 0
    end
    side = 0
    bisection = true
    C = 16 # safety factor for threshold corresponding to 4 iterations = 2^4
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    # calculate k on each bisection step with account for local function properties and symmetry
    for i in 1:maxIter
        local x3, y3
        if bisection
            x3 = (x1 + x2) / 2
            y3 = f(x3) - target  # Function value at midpoint
            ym = (y1 + y2) / 2 # Ordinate of chord at midpoint
            r = 1 - abs(ym / (y2 - y1)) # Symmetry factor
            k = r * r # Deviation factor
            # Check if the function is close enough to linear
            if abs(ym - y3) < k * (abs(y3) + abs(ym))
                bisection = false
                threshold = (x2 - x1) * C
            end
        else
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if x3 <= x1
                x3 = x1
                y3 = y1
            elseif x3 >= x2
                x3 = x2
                y3 = y2
            else
                y3 = f(x3) - target
            end
            threshold /= 2
        end

        if abs(y3) <= eps1 || abs(x2 - x1) <= eps2 # Convergence check
            return x3
        end

        if sign(y1) == sign(y3)
            if side == 1
                m = 1 - y3 / y1
                if m <= 0
                    y2 /= 2
                else
                    y2 *= m
                end
            elseif !bisection
                side = 1
            end
            x1, y1 = x3, y3
        else
            if side == -1
                m = 1 - y3 / y2
                if m <= 0
                    y1 /= 2
                else
                    y1 *= m
                end
            elseif !bisection
                side = -1
            end
            x2, y2 = x3, y3
        end
        if x2 - x1 > threshold # in case AB failed to shrink the interval enough
            bisection = true
            side = 0
        end
    end
    return NaN
end

## --- non bracketing
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
    #Base.depwarn("`secant_method(f, xs)` is deprecated; use `find_zero(f, xs, Secant())` instead.", :secant_method)
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
    #Base.depwarn("`secant(f, a, b)` is deprecated; use `find_zero(f, (a,b), Secant())` instead.", :secant)
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
    #Base.depwarn("`muller(f, x)` is deprecated; use `find_zero(f, x, Roots.Muller())` instead.", :muller)
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

struct TupleWrapper{F,Fp}
    f::F
    fp::Fp
end
(F::TupleWrapper)(x) = begin
    u, v = F.f(x), F.fp(x)
    return (u, u / v)
end

#=
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
=#
newton(f::Tuple, x0; kwargs...) = newton(TupleWrapper(f[1], f[2]), x0; kwargs...)
function newton(f, x0; xatol=nothing, xrtol=nothing, maxevals=100)
    #Base.depwarn("`newton(f, x0)` is deprecated; use `find_zero(f, x0, Roots.Newton())` instead.", :newton)
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
## newton, halley, quadratic_inverse, superhalley, chebyshevlike
"""
    Roots.newton(f, fp, x0; kwargs...)

Implementation of Newton's method: `xᵢ₊₁ =  xᵢ - f(xᵢ)/f'(xᵢ)`.

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- the derivative of `f`.

* `x0::Number` -- initial guess. For Newton's method this may be complex.

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` to be used for the first derivative.

Keyword arguments are passed to `find_zero` using the `Roots.Newton()` method.

See also `Roots.newton((f,fp), x0)` and `Roots.newton(fΔf, x0)` for simpler implementations.

"""
function newton(f, fp, x0; kwargs...)
    #Base.depwarn("`newton(f,fp, x0)` is deprecated; use `find_zero((f,fp), x0, Roots.Newton())` instead.", :newton)
    find_zero((f, fp), x0, Newton(); kwargs...)
end

## --------------------------------------------------
#=
"""
    Roots.halley(f, fp, fpp, x0; kwargs...)

Implementation of Halley's method (cf `?Roots.Halley()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.Halley()` method.

"""
=#
function halley(f, fp, fpp, x0; kwargs...)
    #Base.depwarn("`halley(f,fp,fpp, x0)` is deprecated; use `find_zero((f, fp, fpp), x0, Roots.Halley())` instead.", :halley)
    find_zero((f, fp, fpp), x0, Halley(); kwargs...)
end

#=
"""
    Roots.quadratic_inverse(f, fp, fpp, x0; kwargs...)

Implementation of the quadratic inverse method (cf `?Roots.QuadraticInverse()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `ForwardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.QuadraticInverse()` method.

"""
=#
function quadratic_inverse(f, fp, fpp, x0; kwargs...)
    #Base.depwarn("`quadratic_inverse(f,fp,fpp, x0)` is deprecated; use `find_zero((f, fp, fpp), x0, Roots.QuadraticInverse())` instead.", :quadratic_inverse)
    find_zero((f, fp, fpp), x0, QuadraticInverse(); kwargs...)
end

function superhalley(f, fp, fpp, x0; kwargs...)
    #Base.depwarn("`superhalley(f,fp,fpp, x0)` is deprecated; use `find_zero((f, fp, fpp), x0, Roots.SuperHalley())` instead.", :superhalley)
    find_zero((f, fp, fpp), x0, SuperHalley(); kwargs...)
end

function chebyshev_like(f, fp, fpp, x0; kwargs...)
    #Base.depwarn("`chebyshev_like(f,fp,fpp, x0)` is deprecated; use `find_zero((f, fp, fpp), x0, Roots.ChebyshevLike())` instead.", :chebyshev_like)
    find_zero((f, fp, fpp), x0, ChebyshevLike(); kwargs...)
end
