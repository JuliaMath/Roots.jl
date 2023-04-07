##################################################

# type to throw on successful convergence
mutable struct StateConverged
    x0::Number
end

# type to throw on failure
mutable struct ConvergenceFailed
    reason::AbstractString
end

##################################################
## Helpers for the various methods

_unitless(x) = x / oneunit(x)

# NaN of correct type
nan(::Type{Float16}) = NaN16
nan(::Type{Float32}) = NaN32
nan(::Type{Float64}) = NaN
nan(x::T) where {T<:Number} = NaN * one(T)
nan(x::Type{T}) where {T<:Number} = NaN * one(T)
nan(::Any) = NaN

## issue with approx derivative
isissue(x) = iszero(x) || isnan(x) || isinf(x)

# of (a,fa), (b,fb) choose pair where |f| is smallest
@inline choose_smallest(a, b, fa, fb) = abs(fa) < abs(fb) ? (a, fa) : (b, fb)
@inline sort_smallest(a, b, fa, fb) = abs(fa) < abs(fb) ? (a, b, fa, fb) : (b, a, fb, fa)

# from an interable extract a bracketing interval and promote to floating point values
# used when an interval is desired
_extrema(x::Number) = throw(ArgumentError("Need extrema to return two distinct values"))
function _extrema(x)
    a, b = extrema(x)
    a == b && throw(ArgumentError("Need extrema to return two distinct values"))
    a, b
end
# fix type instability issues of tuples of mixed types
function _extrema(x::Tuple{<:Number,<:Number})
    a, b = x
    a == b && throw(ArgumentError("Need extrema to return two distinct values"))
    extrema(promote(a, b))
end

# used by secant. Get x₀, x₁ for x
function x₀x₁(x::Number)
    x₁ = float(x)
    promote(_default_secant_step(x₁), x₁)
end

function x₀x₁(x)
    x₀, x₁′ = x
    x₁ = first(x₁′)
    promote(float(x₀), float(x₁))
end

## find a default secant step
function _default_secant_step(x1)
    ϵ = eps(one(real(x1)))
    h = cbrt(ϵ)
    dx = h * oneunit(x1) + abs(x1) * h^2 # adjust for if eps(x1) > h
    x0 = x1 + dx
    x0
end

# a bit better than a - fa/f_ab
@inline secant_step(a, b, fa, fb) = a - fa * (b - a) / (fb - fa)

function guarded_secant_step(alpha, beta, falpha, fbeta)
    fp = (fbeta - falpha) / (beta - alpha)
    Δ = fbeta / fp
    ## odd, we get allocations if we define Delta, then beta - Delta
    ## Δ = beta - fbeta * (beta - alpha) / (fbeta - falpha)

    if isissue(Δ)
        Δ = oneunit(alpha) / 1000
    elseif abs(Δ) >= 100 * abs(alpha - beta) # guard runaway
        Δ = sign(Δ) * 100 * min(oneunit(alpha), abs(alpha - beta))
    end

    if isissue(Δ)
        return (alpha + (beta - alpha) * (0.5), true) # midpoint
    else
        return (beta - Δ, false)
    end
end

"""
    steff_step(M, x, fx)

Return first Steffensen step x + fx (with proper units).
May be overridden to provide a guard when fx is too large.

"""
function steff_step(M::Any, x::T, fx::S) where {T,S}
    x + fx * oneunit(T) / oneunit(S)
end

# return vertex of parabola through (a,fa),(b,fb),(c,fc)
# first time through, we have picture of a > b > c; |fa|, |fc| > |fb|, all same sign
function quad_vertex(c, fc, b, fb, a, fa)
    fba = (fb - fa) / (b - a)
    fbc = (fb - fc) / (b - c)

    1 / 2 * ((a + b) - fba / (fbc - fba) * (c - a))
end

## inverse quadratic
function inverse_quadratic_step(a::T, b, c, fa, fb, fc) where {T}
    s = zero(T)
    s += a * fb * fc / (fa - fb) / (fa - fc) # quad step
    s += b * fa * fc / (fb - fa) / (fb - fc)
    s += c * fa * fb / (fc - fa) / (fc - fb)
    s
end

## Different functions for approximating f'(xn)
## return fpxn and whether it is an issue

## use f[a,b] to approximate f'(x)
function _fbracket(a, b, fa, fb)
    num, den = fb - fa, b - a
    iszero(num) && iszero(den) && return Inf, true
    out = num / den
    out, isissue(out)
end

## use f[y,z] - f[x,y] + f[x,z] to approximate
function _fbracket_diff(a, b, c, fa, fb, fc)
    x1, issue = _fbracket(b, c, fb, fc)
    issue && return x1, issue
    x2, issue = _fbracket(a, b, fa, fb)
    issue && return x2, issue
    x3, issue = _fbracket(a, c, fa, fc)
    issue && return x3, issue

    out = x1 - x2 + x3
    out, isissue(out)
end

## use f[a,b] * f[a,c] / f[b,c]
function _fbracket_ratio(a, b, c, fa, fb, fc)
    x1, _ = _fbracket(a, b, fa, fb)
    x2, _ = _fbracket(a, c, fa, fc)
    x3, _ = _fbracket(b, c, fb, fc)
    out = (x1 * x2) / x3
    out, isissue(out)
end

## from https://core.ac.uk/download/pdf/82282744.pdf
## signum based iteration allows one to guarantee a valid starting point
## if N is big enough. (THough this depends on algorithm, a, b and function)
## N here would need to be tuned. But, as is, this may prove useful anyways.
function identify_starting_point(f, a, b, N)
    pts = range(a, stop=b, length=N + 1)
    fxs = f.(pts)
    sfxs = sign.(f.(pts))
    identify_starting_point(a, b, sfxs)
end

function identify_starting_point(a, b, sfxs)
    N = length(sfxs) - 1
    p0 = a + (b - a) / 2
    p1 = p0 + (b - a) / (2N) * sfxs[1] * sum(s for s in sfxs[2:(end - 1)])
    p1
end

## not used
function _unicode_subscript(io, j)
    a = ("⁻", "", "", "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
    for i in string(j)
        print(io, a[Int(i) - 44])
    end
end

unicode_subscript(io, j) = _unicode_subscript.(Ref(io), reverse(digits(j)))
