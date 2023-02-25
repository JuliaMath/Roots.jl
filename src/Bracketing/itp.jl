"""
    Roots.ITP(;[κ₁-0.2, κ₂=2, n₀=1])

Use the [ITP](https://en.wikipedia.org/wiki/ITP_method) bracketing
method.  This method claims it "is the first root-finding algorithm
that achieves the superlinear convergence of the secant method
while retaining the optimal worst-case performance of the bisection
method."

The values `κ1`, `κ₂`, and `n₀` are tuning parameters.

The
[suggested](https://docs.rs/kurbo/0.8.1/kurbo/common/fn.solve_itp.html)
value of `κ₁` is `0.2/(b-a)`, but the default here is `0.2`. The value
of `κ₂` is `2`, and the default value of `n₀` is `1`.

## Note:

Suggested on
[discourse](https://discourse.julialang.org/t/julia-implementation-of-the-interpolate-truncate-project-itp-root-finding-algorithm/77739)
by `@TheLateKronos`, who supplied the original version of the code.

"""
struct ITP{T,S} <: AbstractBracketingMethod
    κ₁::T
    κ₂::S
    n₀::Int
    function ITP(κ₁::T′, κ₂::S, n₀::Int) where {T′,S}
        0 ≤ κ₁ < Inf || throw(ArgumentError("κ₁ must be between 0 and ∞"))
        1 ≤ κ₂ < (3 + √5) / 2 ||
            throw(ArgumentError("κ₂ must be between 1 and 1 plus the golden ratio"))
        0 < n₀ < Inf || throw(ArgumentError("n₀ must be between 0 and ∞"))
        T = float(T′)

        ## κ₂ == 2 || throw(ArgumentError("κ₂ is hardcoded to be 2"))

        new{T,S}(float(κ₁), κ₂, n₀)
    end
end
ITP(; κ₁=0.2, κ₂=2, n₀=1) = ITP(κ₁, κ₂, n₀)

struct ITPState{T,S,R} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    j::Int
    ϵ2n₁₂::R
    d::T
end

function init_state(M::ITP, F, x₀, x₁, fx₀, fx₁)
    if x₀ > x₁
        x₀, x₁, fx₀, fx₁ = x₁, x₀, fx₁, fx₀
    end

    ## we compute this once the options and initial state are known
    ϵ2n₁₂ = zero(float(x₁) / x₁) # ϵ*2^(ceil(Int, log2((b-a)/(2*ϵ))) + n₀)

    # handle interval if fa*fb ≥ 0 (explicit, but also not needed)
    (iszero(fx₀) || iszero(fx₁)) && return ITPState(x₁, x₀, fx₁, fx₀, 0, ϵ2n₁₂, x₁)
    assert_bracket(fx₀, fx₁)

    a, b, fa, fb = x₀, x₁, fx₀, fx₁

    ITPState(b, a, fb, fa, 0, ϵ2n₁₂, a)
end

function update_state(M::ITP, F, o::ITPState{T,S,R}, options, l=NullTracks()) where {T,S,R}
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1
    j, ϵ2n₁₂ = o.j, o.ϵ2n₁₂
    κ₁, κ₂ = M.κ₁, M.κ₂

    if iszero(ϵ2n₁₂)
        # we need the options to set the ϵ⋅2^n₁₂ part of r.
        ϵ = max(options.xabstol, max(abs(a), abs(b)) * options.xreltol)
        ϵ2n₁₂ = ϵ * exp2(ceil(Int, log2((b - a) / (2ϵ))) + M.n₀)
        @set! o.ϵ2n₁₂ = ϵ2n₁₂
    end

    Δ = b - a
    x₁₂ = a + Δ / 2  # middle must be (a+b)/2
    r = ϵ2n₁₂ * exp2(-j) - Δ / 2
    δ′ = κ₁ * Δ^κ₂ # a numeric literal for  κ₂ is faster
    δ = δ′ / oneunit(δ′)
    # δ = κ₁ * Δ^2
    xᵣ = (b * fa - a * fb) / (fa - fb)

    σ = sign(x₁₂ - xᵣ)
    xₜ = δ ≤ abs(x₁₂ - xᵣ) / oneunit(xᵣ) ? xᵣ + σ * δ * oneunit(xᵣ) : x₁₂

    c::T = xᵢₜₚ = abs(xₜ - x₁₂) ≤ r ? xₜ : x₁₂ - σ * r

    if !(a < c < b)
        nextfloat(a) ≥ b &&
            log_message(l, "Algorithm stopped narrowing bracketing interval")
        return (o, true)
    end

    fc::S = F(c)
    incfn(l)

    if sign(fa) * sign(fc) < 0
        b, fb = c, fc
    else
        a, fa = c, fc
    end

    o = _set(o, (b, fb), (a, fa))
    @set! o.j = j + 1

    return o, false
end
