"""
    Roots.Ridders()

Implements [Ridders'](https://en.wikipedia.org/wiki/Ridders%27_method) method.
This bracketing method finds the midpoint, `x₁`; then interpolates an exponential; then uses false position with the interpolated value to find `c`. If `c` and `x₁` form a bracket is used, otherwise the subinterval `[a,c]` or `[c,b]` is used.

Example:

```jldoctest
julia> using Roots

julia> find_zero(x -> exp(x) - x^4, (5, 15), Roots.Ridders()) ≈ 8.61316945644
true

julia> find_zero(x -> x*exp(x) - 10, (-100, 100), Roots.Ridders()) ≈ 1.74552800274
true

julia> find_zero(x -> tan(x)^tan(x) - 1e3, (0, 1.5), Roots.Ridders()) ≈ 1.3547104419
true
```

[Ridders](https://cs.fit.edu/~dmitra/SciComp/Resources/RidderMethod.pdf) showed the error satisfies `eₙ₊₁ ≈ 1/2 eₙeₙ₋₁eₙ₋₂ ⋅ (g^2-2fh)/f` for
`f=F', g=F''/2, h=F'''/6`, suggesting converence at rate `≈ 1.839...`. It uses two function evaluations per step, so
 its order of convergence is `≈ 1.225...`.
"""
struct Ridders <: AbstractBracketingMethod end

function update_state(
    M::Ridders,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1

    x₁ = a + (b - a) / 2
    fx₁ = F(x₁)
    incfn(l)

    c::T = x₁ + (x₁ - a) * sign(fa) * fx₁ / sqrt(fx₁^2 - fa * fb)
    fc::S = F(c)
    incfn(l)

    if !(a < c < b)
        nextfloat(a) ≥ b &&
            log_message(l, "Algorithm stopped narrowing bracketing interval")
        return (o, true)
    end

    # choose bracketing interval from [x₁, c], [c, x₁], [a,c], [c,b]
    if sign(fx₁) * sign(fc) < 0
        a, b, fa, fc = x₁ < c ? (x₁, c, fx₁, fc) : (c, x₁, fc, fx₁)
    elseif sign(fa) * sign(fc) < 0
        b, fb = c, fc
    else
        a, fa = c, fc
    end

    o = _set(o, (b, fb), (a, fa))

    return o, false
end
