
### Options

abstract type AbstractUnivariateZeroOptions end

struct UnivariateZeroOptions{Q,R,S,T} <: AbstractUnivariateZeroOptions
    xabstol::Q
    xreltol::R
    abstol::S
    reltol::T
    maxiters::Int
    strict::Bool
    exact::Symbol
end


function init_options(M::AbstractUnivariateZeroMethod,
                      state::AbstractUnivariateZeroState{T,S};
                      kwargs...) where {T,S}
    d = kwargs
    defs = default_tolerances(M, state)
    δₐ = get(d, :xatol, get(d, :xabstol, defs[1]))
    δᵣ = get(d, :xrtol, get(d, :xreltol, defs[2]))
    ϵₐ = get(d, :atol, get(d, :abstol, defs[3]))
    ϵᵣ = get(d, :rtol, get(d, :reltol, defs[4]))
    maxiters = get(d, :maxiters, get(d, :maxevals, get(d, :maxsteps, defs[5])))
    strict = get(d, :strict, defs[6])

    exact = :nothing
    if iszero(δₐ) && iszero(δᵣ)
        exact = :xexact
    elseif iszero(ϵₐ) && iszero(ϵᵣ)
        exact = :fexact
    end
    return UnivariateZeroOptions(δₐ, δᵣ, ϵₐ, ϵᵣ, maxiters, strict, exact)
end

## --------------------------------------------------

"""
    default_tolerances(M::AbstractUnivariateZeroMethod, state::AbstractUnivariateZeroState{T,S})

The default tolerances for most methods are `xatol=eps(T)`,
`xrtol=eps(T)`, `atol=4eps(S)`, and `rtol=4eps(S)`, with the proper
units (absolute tolerances have the units of `x` and `f(x)`; relative
tolerances are unitless). For `Complex{T}` values, `T` is used.

The number of iterations is limited by `maxiters=40`.

"""
function default_tolerances(
    ::AbstractUnivariateZeroMethod,
    ::AbstractUnivariateZeroState{T, S}
) where {T,S}
    xatol = eps(real(T)) * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4 * eps(real(float(S))) * one(real(S))
    maxiters = 40
    strict = false
    (xatol, xrtol, atol, rtol, maxiters, strict)
end

## --------------------------------------------------

# ## Assess convergence

## test f == 0 not f ≈ 0
function is_exact_zero_f(
    ::AbstractNonBracketingMethod,
    state::AbstractUnivariateZeroState,
    options,
)
    fb = state.fxn1
    iszero(fb)
end

function is_exact_zero_f(
    ::AbstractBracketingMethod,
    state::AbstractUnivariateZeroState,
    options,
)
    fa, fb = state.fxn0, state.fxn1
    iszero(fa) || iszero(fb)
end

## test f ≈ 0 not f == 0
function is_approx_zero_f(
    ::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState,
    options::O,
) where {O<:AbstractUnivariateZeroOptions}
    ab, afb = abs(state.xn1), abs(state.fxn1)
    ϵₐ, ϵᵣ = options.abstol, options.reltol
    Δ = max(_unitless(ϵₐ), _unitless(ab) * ϵᵣ)
    afb ≤ Δ * oneunit(afb)
end

## test f ≈ 0 not f == 0
function is_approx_zero_f(
    ::AbstractBracketingMethod,
    state::AbstractUnivariateZeroState,
    options::O,
) where {O<:AbstractUnivariateZeroOptions}
    ab₁, afb₁ = abs(state.xn1), abs(state.fxn1)
    ab₀, afb₀ = abs(state.xn0), abs(state.fxn0)
    ϵₐ, ϵᵣ = options.abstol, options.reltol
    u, fu = afb₀ < afb₁ ? (ab₀, afb₀) : (ab₁, afb₁)
    Δ = max(_unitless(ϵₐ), _unitless(u) * ϵᵣ)
    fu ≤ Δ * oneunit(fu)
end

function is_approx_zero_f(
    ::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState,
    options::O,
    relaxed::Any,
) where {O<:AbstractUnivariateZeroOptions}
    ab, afb = abs(state.xn1), abs(state.fxn1)
    ϵₐ, ϵᵣ = options.abstol, options.reltol
    Δ = max(_unitless(ϵₐ), _unitless(ab) * ϵᵣ)
    Δ = cbrt(abs(_unitless(Δ))) * oneunit(afb) # relax test
    afb <= Δ
end

## --------------------------------------------------

# testing xₙ₊₁ - xₙ ≈ 0
function iszero_Δx(
    ::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState,
    options
)
    a, b, fa, fb = state.xn0, state.xn1, state.fxn0, state.fxn1
    δₐ, δᵣ = options.xabstol, options.xreltol
    isapprox(a, b, atol=δₐ, rtol=δᵣ)
end

# this is for Bisection where we can go adjacent floating point values
function iszero_Δx_xexact(
    ::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState,
    options,
)
    a, b = state.xn0, state.xn1
    if b < a
        a, b = b, a
    end

    nextfloat(float(a)) == float(b)
end

# this is Alefeld, Potra, Shi criteria for termination
function iszero_Δx_fexact(
    ::AbstractBracketingMethod,
    state::AbstractUnivariateZeroState,
    options
)
    a, b, fa, fb = state.xn0, state.xn1, state.fxn0, state.fxn1
    u, fu = choose_smallest(a, b, fa, fb)
    δₐ, δᵣ = options.xabstol, options.xreltol
    δₓ = max(δₐ, 2 * (abs(u) * δᵣ)) # needs non-zero δₐ to stop near 0
    abs(b - a) ≤ δₓ
end

function iszero_Δx_exact(M::AbstractBracketingMethod,
                         state::AbstractUnivariateZeroState{T,S},
                         options) where {T, S}
    if T ∈ (Float16, Float32, Float64)
        iszero_Δx_xexact(M, state, options)
    else
        iszero_Δx_fexact(M, state, options)
    end
end

# test when fconverged to ensure not runaway
function is_small_Δx(
    M::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState,
    options
)
    δ = _unitless(abs(state.xn1 - state.xn0))
    δₐ, δᵣ = options.xabstol, options.xreltol
    Δₓ = max(_unitless(δₐ), _unitless(abs(state.xn1)) * δᵣ)
    Δₓ = sqrt(sqrt(sqrt((abs(_unitless(Δₓ)))))) # faster than x^(1/8)
    return δ ≤ Δₓ
end

isnan_x(M::AbstractBracketingMethod, state) = isnan(state.xn1) || isnan(state.xn0)
isnan_x(M::AbstractNonBracketingMethod, state) = isnan(state.xn1)

isinf_x(M::AbstractBracketingMethod, state) = isinf(state.xn1) || isinf(state.xn0)
isinf_x(M::AbstractNonBracketingMethod, state) = isinf(state.xn1)

isnan_f(M::AbstractBracketingMethod, state) = isnan(state.fxn1) || isnan(state.fxn0)
isnan_f(M::AbstractNonBracketingMethod, state) = isnan(state.fxn1)

isinf_f(M::AbstractBracketingMethod, state) = isinf(state.fxn1) || isinf(state.fxn0)
isinf_f(M::AbstractNonBracketingMethod, state) = isinf(state.fxn1)

## --------------------------------------------------

"""
    Roots.assess_convergence(method, state, options)

Assess if algorithm has converged.

Return a convergence flag and a Boolean indicating if algorithm has terminated (converged or not converged)

If algorithm hasn't converged this returns `(:not_converged, false)`.

If algorithm has stopped or converged, return flag and `true`. Flags are:

* `:x_converged` if `xn1 ≈ xn`, typically with non-zero tolerances specified.

* `:f_converged` if  `|f(xn1)| < max(atol, |xn1|*rtol)`

* `:nan` or `:inf` if fxn1 is `NaN` or an infinity.

* `:not_converged` if algorithm should continue

Does not check number of steps taken nor number of function evaluations.

In `decide_convergence`, stopped values (and `:x_converged` when `strict=false`) are checked for convergence with a relaxed tolerance.


"""
function assess_convergence(M::Any, state::AbstractUnivariateZeroState, options)
    _assess_convergence(M, state, options)
end

function _assess_convergence(M, state, options)
    # return convergence_flag, boolean
    is_exact_zero_f(M, state, options) && return (:exact_zero, true)
    isnan_x(M, state) && return (:nan, true)
    isinf_x(M, state) && return (:inf, true)
    isnan_f(M, state) && return (:nan, true)
    isinf_f(M, state) && return (:inf, true)
    is_approx_zero_f(M, state, options) && return (:f_converged, true)
    iszero_Δx(M, state, options) && return (:x_converged, true)
    return (:not_converged, false)
end

function assess_convergence(
    M::AbstractBracketingMethod,
    state::AbstractUnivariateZeroState,
    options
)
    if options.exact == :xexact || options.exact == :exact
        # return convergence_flag, boolean
        # no check if f == ∞
        is_exact_zero_f(M, state, options)  && return (:exact_zero, true)
        isnan_f(M, state)                   && return (:nan, true)
        is_approx_zero_f(M, state, options) && return (:f_converged, true)
        iszero_Δx_exact(M, state, options)  && return (:x_converged, true)
    elseif options.exact == :fexact
        is_exact_zero_f(M, state, options)  && return (:exact_zero, true)
        isnan_f(M, state)                   && return (:nan, true)
        iszero_Δx_fexact(M, state, options) && return (:x_converged, true)
    else
        return _assess_convergence(M, state, options)
    end

    return (:not_converged, false)
end


#=
"""
    decide_convergence(M,F,state,options, convergence_flag)

When the algorithm terminates, this function decides the stopped value or returns NaN
"""
=#
function decide_convergence(
    M::AbstractNonBracketingMethod,
    F,
    state::AbstractUnivariateZeroState{T,S},
    options,
    val,
) where {T,S}
    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    val ∈ (:f_converged, :exact_zero, :converged) && return xn1

    ## XXX this could be problematic
    if val == :nan
        # return if Δx small
        Δₓ = abs(xn1 - xn0)
        δₐ, δᵣ = options.xabstol, options.xreltol
        u = min(abs(xn0), abs(xn1))
        δₓ = max(δₐ, 2 * (abs(u) * δᵣ)) # needs non-zero δₐ to stop near 0
        Δₓ ≤ δₓ && return xn1

        # or if abs(fxn0) small
        ϵₐ, ϵᵣ = options.abstol, options.reltol
        Δ = max(_unitless(ϵₐ), _unitless(xn0) * ϵᵣ)
        abs(state.fxn0) ≤ Δ * oneunit(state.fxn0) && return xn0

        # else
        return nan(T) * xn1
    end
    val == :inf_nan && return xn1

    ## stopping is a heuristic, x_converged can mask issues
    ## if strict=true or the tolerance for f is 0 this will return xn1 if x_converged
    ## if strict == false, this will also check f(xn) ~ - with a relaxed
    ## tolerance
    fexact = options.exact == :fexact
    if options.strict || fexact
        val == :x_converged && return xn1
        !fexact && is_approx_zero_f(M, state, options) && return xn1
        fexact  && is_exact_zero_f(M, state, options) && return xn1
    else
        if val == :x_converged
            # The :xexact case isn't always spelled out in the type, so
            # we replicate a bit here
            δ, ϵ = options.abstol, options.reltol
            if (iszero(δ) && iszero(ϵ)) && (options.exact != :xexact)
                error("Hah")
            end

            iszero(δ) && iszero(ϵ) && return xn1
            is_approx_zero_f(M, state, options, true) && return xn1
        elseif val == :not_converged
            # this is the case where runaway can happen
            ## XXX Need a good heuristic to catch that
            is_approx_zero_f(M, state, options, :relaxed) &&
                is_small_Δx(M, state, options) &&
                return xn1
        end
    end

    nan(T) * xn1
end

# assumes stopped = :x_converged
function decide_convergence(
    M::AbstractBracketingMethod,
    F,
    state::AbstractUnivariateZeroState,
    options,
    val,
)
    a, b = state.xn0, state.xn1
    fa, fb = state.fxn0, state.fxn1

    iszero(fa) && return a
    iszero(fb) && return b
    isnan(fa) && return a
    isnan(fb) && return b

    # get as close as possible with one extra function call
    # when exact closeness is possible
    a₊₊ = nextfloat(nextfloat(float(a)))
    if b == a₊₊
        c = nextfloat(float(a))
        fc = first(F(c))
        m = minimum(abs, (fa, fb, fc))
        abs(fc) == m && return c
        abs(fa) == m && return a
        return b
    end

    return (abs(fa) < abs(fb)) ? a : b
end
