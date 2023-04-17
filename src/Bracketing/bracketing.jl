### Bracketing method defaults

function init_state(M::AbstractBracketingMethod, F::Callable_Function, x)
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function init_state(M::AbstractBracketingMethod, F, x₀, x₁, fx₀, fx₁)
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)
    assert_bracket(fx₀, fx₁)
    a, b, fa, fb = (x₀ < x₁) ? (x₀, x₁, fx₀, fx₁) : (x₁, x₀, fx₁, fx₀)
    UnivariateZeroState(b, a, fb, fa)
end

Base.last(state::AbstractUnivariateZeroState, M::AbstractBracketingMethod) =
    state.xn0 < state.xn1 ? (state.xn0, state.xn1) : (state.xn1, state.xn0)

fn_argout(::AbstractBracketingMethod) = 1
initial_fncalls(::AbstractBracketingMethod) = 2

## tracks for bisection, different from secant, we show bracketing interval
## No init here; for Bisection() [a₀, b₀] is just lost.
function log_step(l::Tracks, M::AbstractBracketingMethod, state; init::Bool=false)
    a, b = state.xn0, state.xn1
    push!(l.abₛ, a < b ? (a, b) : (b, a))
    !init && log_iteration(l, 1)
    nothing
end

# use xatol, xrtol only, but give some breathing room over the strict ones and cap number of steps
function default_tolerances(::AbstractBracketingMethod, ::Type{T}, ::Type{S}) where {T,S}
    xatol = eps(real(T))^3 * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 0 * oneunit(real(S))
    rtol = 0 * one(real(S))
    maxevals = 60
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, strict)
end

## --------------------------------------------------

const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""

## utils
@inline isbracket(fa, fb) = sign(fa) * sign(fb) < 0
assert_bracket(fx0, fx1) = isbracket(fx0, fx1) || throw(ArgumentError(bracketing_error))

## helper function: floating point, sorted, finite
function adjust_bracket(x0)
    u, v = map(float, _extrema(x0))
    if u > v
        u, v = v, u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end
