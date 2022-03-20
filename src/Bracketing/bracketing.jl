###

const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""

## utils
@inline isbracket(fa, fb) = sign(fa) * sign(fb) < 0
assert_bracket(fx0, fx1) = isbracket(fx0, fx1) || throw(ArgumentError(bracketing_error))


## tracks for bisection, different, we show bracketing interval
## No init here; for Bisection() [a₀, b₀] is just lost.
function log_step(l::Tracks, M::AbstractBracketing, state; init::Bool=false)
    a, b = state.xn0, state.xn1
    push!(l.abₛ, a < b ? (a,b) : (b,a))
    !init && log_iteration(l, 1)
    nothing
end


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

function init_state(M::AbstractBracketing, F::Callable_Function, x)
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = F(x₀), F(x₁)
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function init_state(::AbstractBracketing, F, x₀, x₁, fx₀, fx₁; m=_middle(x₀, x₁), fm=F(m))
    if x₀ > x₁
        x₀, x₁, fx₀, fx₁ = x₁, x₀, fx₁, fx₀
    end

    # handle interval if fa*fb ≥ 0 (explicit, but also not needed)
    (iszero(fx₀) || iszero(fx₁)) && return UnivariateZeroState(x₁, x₀, fx₁, fx₀)
    assert_bracket(fx₀, fx₁)

    if sign(fm) * fx₀ < 0
        a, b, fa, fb = x₀, m, fx₀, fm
    else
        a, b, fa, fb = m, x₁, fm, fx₁
    end

    sign(a) * sign(b) < 0 && throw(ArgumentError("_middle error"))

    UnivariateZeroState(b, a, fb, fa)
end


initial_fncalls(::Roots.AbstractBracketing) = 2
fn_argout(::AbstractBracketing) = 1



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

    M = maximum(abs, (a, b))
    δₓ = maximum(
        promote(options.xabstol, M * options.xreltol, sign(options.xreltol) * eps(M)),
    )

    if abs(b - a) <= 2δₓ
        return (:x_converged, true)
    end

    # check f
    u, fu = choose_smallest(a, b, fa, fb)
    δ = maximum(promote(options.abstol, M * options.reltol * (oneunit(fu) / oneunit(u))))
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

    if abs(fa) < abs(fb)
        return a
    else
        return b
    end
end


## convergence is much different here
function check_zero(::AbstractBracketing, state, c, fc)
    isnan(c) && return true
    isinf(c) && return true
    iszero(fc) && return true
    return false
end
