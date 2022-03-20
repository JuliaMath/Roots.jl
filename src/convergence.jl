### Options
struct UnivariateZeroOptions{Q,R,S,T}
    xabstol::Q
    xreltol::R
    abstol::S
    reltol::T
    maxevals::Int
    maxfnevals::Int
    strict::Bool
end


init_options(
    M::AbstractUnivariateZeroMethod,
    state::AbstractUnivariateZeroState{T,S};
    kwargs...,
) where {T,S} = init_options(M, T, S; kwargs...)

function init_options(M, T=Float64, S=Float64; kwargs...)
    d = kwargs

    defs = default_tolerances(M, T, S)
    options = UnivariateZeroOptions(
        get(d, :xatol, get(d, :xabstol, defs[1])),
        get(d, :xrtol, get(d, :xreltol, defs[2])),
        get(d, :atol, get(d, :abstol, defs[3])),
        get(d, :rtol, get(d, :reltol, defs[4])),
        get(d, :maxevals, get(d, :maxsteps, defs[5])),
        get(d, :maxfnevals, defs[6]),
        get(d, :strict, defs[7]),
    )
    if haskey(d, :maxfnevals)
        @warn(
            "maxfnevals is ignored. See the test for an example to implement this featrue"
        )
    end
    options
end

# # reset options to default values
@deprecate init_options!(options, M) init_options(M)

## --------------------------------------------------

"""
    default_tolerances(M::AbstractUnivariateZeroMethod, [T], [S])

The default tolerances for most methods are `xatol=eps(T)`,
`xrtol=eps(T)`, `atol=4eps(S)`, and `rtol=4eps(S)`, with the proper
units (absolute tolerances have the units of `x` and `f(x)`; relative
tolerances are unitless). For `Complex{T}` values, `T` is used.

The number of iterations is limited by `maxevals=40`.

"""
default_tolerances(M::AbstractUnivariateZeroMethod) =
    default_tolerances(M, Float64, Float64)
function default_tolerances(
    ::AbstractUnivariateZeroMethod,
    ::Type{T},
    ::Type{S},
) where {T,S}
    xatol = eps(real(T)) * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4 * eps(real(float(S))) * one(real(S))
    maxevals = 40
    maxfnevals = typemax(Int)
    strict = false
    (xatol, xrtol, atol, rtol, maxevals, maxfnevals, strict)
end

## --------------------------------------------------

## Assess convergence
@inline function _is_f_approx_0(fa, a, atol, rtol, relaxed::Any)
    aa, afa = abs(a), abs(fa)
    tol = max(_unitless(atol), _unitless(aa) * rtol)
    tol = cbrt(abs(_unitless(tol)))  # relax test
    afa <= tol * oneunit(afa)
end
@inline function _is_f_approx_0(fa, a, atol, rtol)
    aa, afa = abs(a), abs(fa)
    tol = max(_unitless(atol), _unitless(aa) * rtol)
    afa <= tol * oneunit(afa)
end

"""
    Roots.assess_convergence(method, state, options)

Assess if algorithm has converged.

Return a convergence flag and a Boolean indicating if algorithm has terminated (converged or not converged)

If algrithm hasn't converged returns `(:not_converged, false)`.

If algorithm has stopped or converged, return flag and `true`. Flags are:

* `:x_converged` if `abs(xn1 - xn0) < max(xatol, max(abs(xn1), abs(xn0)) * xrtol)`

* `:f_converged` if  `|f(xn1)| < max(atol, |xn1|*rtol)`

* `:nan`, `:inf` if xn1 or fxn1 is `NaN` or an infinity

* `:not_converged` if algorithm should continue

Does not check number of steps taken nor number of function evaluations.

In `decide_convergence`, stopped values (and `:x_converged` when `strict=false`) are checked for convergence with a relaxed tolerance.


"""
function assess_convergence(::Any, state::AbstractUnivariateZeroState, options)
    # return convergence_flag, boolean

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    if isnan(xn1) || isnan(fxn1)
        return (:nan, true)
    end

    if isinf(xn1) || isinf(fxn1)
        return (:inf, true)
    end

    # f(xstar) ≈ xstar * f'(xstar)*eps(), so we pass in lambda
    if _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol)
        return (:f_converged, true)
    end

    # stop when xn1 ~ xn.
    # in find_zeros there is a check that f could be a zero with a relaxed tolerance
    if isapprox(xn1, xn0, atol=options.xabstol, rtol=options.xreltol)
        return (:x_converged, true)
    end

    return (:not_converged, false)
end

# state has stopped, this identifies if it has converged
"""
    decice_convergence(M,F,state,options, convergence_flag)

When the algorithm terminates, this function decides the stopped value or returns NaN
"""
function decide_convergence(
    ::AbstractUnivariateZeroMethod,
    F,
    state::AbstractUnivariateZeroState,
    options,
    val,
)
    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    val ∈ (:f_converged, :exact_zero, :converged) && return xn1

    ## stopping is a heuristic, x_converged can mask issues
    ## if strict=true or the tolerance for f is 0 this will return xn1 if x_converged
    ## if strict == false, this will also check f(xn) ~ - with a relaxed
    ## tolerance
    if options.strict || (iszero(options.abstol) && iszero(options.reltol))
        val == :x_converged && return xn1
        _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol) && return xn1
    else
        if val ∈ (:x_converged, :not_converged)
            # |fxn| <= relaxed

            _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol, true) && return xn1
        end
    end

    # if (state.stopped || state.x_converged) && !(state.f_converged)
    #     ## stopped is a heuristic, x_converged can mask issues
    #     ## if strict == false, this will also check f(xn) ~ - with a relaxed
    #     ## tolerance

    #     ## are we at a crossing values?
    #     ## seems worth a check for 2 fn evals.
    #     if T <: Real && S <: Real
    #         for u in (prevfloat(xn1), nextfloat(xn1))
    #             fu = first(F(u))
    #             incfn(state)
    #             if iszero(fu) || _unitless(fu * fxn1) < 0
    #                 state.message *= "Change of sign at xn identified. "
    #                 state.f_converged = true
    #             end
    #         end
    #     end

    #     δ = maximum(_unitless, (options.abstol, options.reltol))
    #     if options.strict || iszero(δ)
    #         if state.x_converged
    #             state.f_converged = true
    #         else
    #             state.convergence_failed = true
    #         end

    #     else
    #         xstar, fxstar = state.xn1, state.fxn1
    #         if _is_f_approx_0(fxstar, xstar, options.abstol, options.reltol, :relaxed)
    #             state.xstar, state.fxstar = xstar, fxstar
    #             msg = "Algorithm stopped early, but |f(xn)| < ϵ^(1/3), where ϵ depends on xn, rtol, and atol. "
    #             state.message = state.message == "" ? msg : state.message * "\n\t" * msg
    #             state.f_converged = true
    #         else
    #             state.convergence_failed = true
    #         end
    #     end
    # end

    # if !state.f_converged
    #     state.xstar, state.fxstar = NaN*xn1, NaN*fxn1
    # end

    NaN * xn1
end
