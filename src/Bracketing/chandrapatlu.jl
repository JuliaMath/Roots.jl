"""
    Roots.Chandrapatla()

Use [Chandrapatla's
algorithm](https://doi.org/10.1016/S0965-9978(96)00051-8)
(cf. [Scherer](https://www.google.com/books/edition/Computational_Physics/cC-8BAAAQBAJ?hl=en&gbpv=1&pg=PA95&printsec=frontcover))
to solve ``f(x) = 0``.

Chandrapatla's algorithm chooses between an inverse quadratic step or a bisection step based on a computed inequality.


"""
struct Chandrapatla <: AbstractBracketingMethod end

struct ChandrapatlaState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    c::T  # keep xₙ₋₂ around for quadratic step
    fxn1::S
    fxn0::S
    fc::S
end

# a = most recent, b prior
function init_state(::Chandrapatla, F, x₀, x₁, fx₀, fx₁)
    a, b, fa, fb = x₁, x₀, fx₁, fx₀
    c, fc = a, fa
    ChandrapatlaState(a, b, c, fa, fb, fc)
end

function update_state(
    ::Chandrapatla,
    F,
    o::ChandrapatlaState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    a, b, c = o.xn1, o.xn0, o.c
    fa, fb, fc = o.fxn1, o.fxn0, o.fc

    # encoding: a = xₙ, b=xₙ₋₁, c= xₙ₋₂

    ξ = (a - b) / (c - b)
    ϕ = (fa - fb) / (fc - fb)
    ϕ² = ϕ^2
    Δ = (ϕ² < ξ) && (1 - 2ϕ + ϕ² < 1 - ξ) # Chandrapatla's inequality to determine next step

    xₜ::T = Δ ? inverse_quadratic_step(a, b, c, fa, fb, fc) : _middle(a, b)

    fₜ::S = F(xₜ)
    incfn(l)

    if sign(fₜ) * sign(fa) < 0
        a, b, c = xₜ, a, b
        fa, fb, fc = fₜ, fa, fb
    else
        a, c = xₜ, a
        fa, fc = fₜ, fa
    end

    o = _set(o, (a, fa), (b, fb)) # a is xₙ, b is xₙ₋₁
    @set! o.c = c
    @set! o.fc = fc

    return (o, false)
end
