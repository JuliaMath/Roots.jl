"""
    Roots.Order1B()
    Roots.King()

A superlinear (order `1.6...`) modification of the secant method for multiple roots.
Presented in A SECANT METHOD FOR MULTIPLE ROOTS, by RICHARD F. KING, BIT 17 (1977), 321-328

The basic idea is similar to Schroder's method: apply the secant method
to  `f/f'`. However, this uses `f' ~ fp = (fx - f(x-fx))/fx` (a Steffensen step). In
this implementation, `Order1B`, when `fx` is too big, a single secant step of `f`
is used.

The *asymptotic* error, `eᵢ = xᵢ - α`, is given by
`eᵢ₊₂ = 1/2⋅G''/G'⋅ eᵢ⋅eᵢ₊₁ + (1/6⋅G'''/G' - (1/2⋅G''/G'))^2⋅eᵢ⋅eᵢ₊₁⋅(eᵢ+eᵢ₊₁)`.

"""
struct King <: AbstractSecantMethod end
struct Order1B <: AbstractSecantMethod end

struct KingState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    G0::S
end

function init_state(::Union{King,Order1B}, F, x₀, x₁, fx₀, fx₁)
    fₛ₀ = F(x₀ - fx₀ * oneunit(x₀) / oneunit(fx₀))
    G₀ = -fx₀^2 / (fₛ₀ - fx₀)
    KingState(x₁, x₀, fx₁, fx₀, G₀)
end
initial_fncalls(::Union{King,Order1B}) = 3

function update_state(::Order1B, F, o::KingState, options, l=NullTracks())
    if do_guarded_step(Order1B(), o)
        state, flag = update_state(Order1(), F, o, options, l)

        x0, fx0 = state.xn0, state.fxn0 # clunky! Need to update G₀ after Order1() step
        fₛ = F(x0 - fx0 * oneunit(x0) / oneunit(fx0))
        incfn(l)
        G₀ = -fx0^2 / (fₛ - fx0)
        @set! state.G0 = G₀

        return (state, flag)
    else
        update_state(King(), F, o, options, l)
    end
end

function update_state(::King, F, o::KingState{T,S}, options, l=NullTracks()) where {T,S}
    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1
    G₀ = o.G0
    fₛ₁ = F(x1 - fx1 * oneunit(x1) / oneunit(fx1))
    incfn(l)
    G₁ = -fx1^2 / (fₛ₁ - fx1)

    m = (x1 - x0) / (G₁ - G₀) # approximate value of `m`, the multiplicity

    if abs(m) <= 1e-2 * oneunit(m)
        log_message(l, "Estimate for multiplicity has issues. ")
        return (o, true)
    end

    Δ = G₁ * (x1 - x0) / (G₁ - G₀)

    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1::T = x1, x1 - Δ
    fx0, fx1::S = fx1, F(x1)
    incfn(l)

    o = _set(o, (x1, fx1), (x0, fx0))
    @set! o.G0 = G₁

    return o, false
end
