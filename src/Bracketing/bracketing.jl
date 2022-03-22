### Bracketing method defaults

function init_state(M::AbstractBracketing, F::Callable_Function, x)
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = F(x₀), F(x₁)
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function init_state(M::AbstractBracketing, F, x₀, x₁, fx₀, fx₁)
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)
    assert_bracket(fx₀, fx₁)
    a, b, fa, fb = (x₀ < x₁) ? (x₀, x₁, fx₀, fx₁) : (x₁, x₀, fx₁, fx₀)
    UnivariateZeroState(b, a, fb, fa)
end

Base.last(state::AbstractUnivariateZeroState, M::AbstractBracketing) =
    state.xn0 < state.xn1 ? (state.xn0, state.xn1) : (state.xn1, state.xn0)

fn_argout(::AbstractBracketing) = 1
initial_fncalls(::AbstractBracketing) = 2

## tracks for bisection, different from secant, we show bracketing interval
## No init here; for Bisection() [a₀, b₀] is just lost.
function log_step(l::Tracks, M::AbstractBracketing, state; init::Bool=false)
    a, b = state.xn0, state.xn1
    push!(l.abₛ, a < b ? (a,b) : (b,a))
    !init && log_iteration(l, 1)
    nothing
end

# use xatol, xrtol only, but give some breathing room over the strict ones and cap number of steps
function default_tolerances(::AbstractBracketing, ::Type{T}, ::Type{S}) where {T,S}
    xatol = eps(real(T))^3 * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 0 * oneunit(real(S))
    rtol = 0 * one(real(S))
    maxevals = 60
    maxfnevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end



function assess_convergence(
    method::AbstractBracketing,
    state::AbstractUnivariateZeroState,
    options,
)
    a, b, fa, fb = state.xn0, state.xn1, state.fxn0, state.fxn1

    if isnan(a) || isnan(b)
        return (:nan, true)
    end

    if isnan(fa) || isnan(fb)
        return (:nan, true)
    end

    # check |b-a| ≤ 2 |u| ϵ + ϵₐ where u ∈ {a,b} is chosen the smaller of |f(a)|, |f(b)|
    u, fu = choose_smallest(a, b, fa, fb)
    δₓ = max(options.xabstol, 2 * abs(u) * options.xreltol) # needs non-zero xabstol to stop near 0
    abs(b-a) ≤ δₓ && return (:x_converged, true)

    # check f (typically not used!)
    δ = max(options.abstol, (u / oneunit(u)) * (options.reltol * oneunit(fu)))

    if abs(fu) <= δ
        iszero(fu) && return (:exact_zero, true)
        return (:f_converged, true)
    end

    return (:not_converged, false)
end


# assumes stopped = :x_converged
function decide_convergence(
    ::AbstractBracketing,
    F,
    state::AbstractUnivariateZeroState,
    options,
    val,
)
    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    isnan(fa) && return a
    isnan(fb) && return b


    abs(fa) < abs(fb) ? a : b
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
