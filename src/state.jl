### States
struct UnivariateZeroState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
end

# simple helper to set main properties of a state object
function _set(state, xf1)
    x, fx = xf1
    @set! state.xn1 = x
    @set! state.fxn1 = fx

    state
end

function _set(state, xf1, xf0)
    x, fx = xf1
    @set! state.xn1 = x
    @set! state.fxn1 = fx

    x, fx = xf0
    @set! state.xn0 = x
    @set! state.fxn0 = fx

    state
end

# init_state(M, F, x; kwargs...)
# init_state(M, F x₀,x₁,fx₀,fx₁; kwargs...)
# init_state(M, state, F)
#
# A state holds at a minimum:
#
# * the values xₙ₋₁, xₙ and f(xₙ₋₁), f(xₙ) along with
# * some method-specific values
#
#
# A state is initialized with `init_state(M, F, x)` which sets up xₙ₋₁, xₙ, f(xₙ₋₁), f(xₙ)
# which then calls `init_state(M, F, xₙ₋₁, xₙ, f(xₙ₋₁), f(xₙ))` to finish the initialization
# to change to a new state use `init_state(M, state, F)`

# basic idea to convert from N to M:
# Fₘ = some state
# stateₘ = init_state(M, stateₙ, Fₘ)
function init_state(M::AbstractUnivariateZeroMethod, state::AbstractUnivariateZeroState, F)
    init_state(M, F, state.xn0, state.xn1, state.fxn0, state.fxn1)
end

# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractUnivariateZeroMethod, F, x)
    error("no default method")
end

# initialize from xs, fxs
function init_state(::AbstractUnivariateZeroMethod, F, x₀, x₁, fx₀, fx₁)
    error("no default method")
end

Base.last(state::AbstractUnivariateZeroState, M::AbstractNonBracketingMethod) = state.xn1
