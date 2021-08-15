using Setfield


struct SSecant <: AbstractSecant end
export SSecant
fn_argout(::SSecant) = 1

struct SNewton <: AbstractSecant end
export SNewton
fn_argout(::SNewton) = 2

struct SBisection <: AbstractBisection end
export SBisection
fn_argout(::SBisection) = 1


abstract type StaticUnivariateZeroState  <:  AbstractUnivariateZeroState end
struct SUnivariateZeroState{T,S} <: StaticUnivariateZeroState where {T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
end

init_options(M::AbstractUnivariateZeroMethod,
             state::StaticUnivariateZeroState;
             kwargs...
             )  = init_options(M, TS(state)...; kwargs...)


# init_state(M, F, state) -- convert
# basic idea to convert from N to M:
# Fₘ = copy(M, fₙ)
# stateₘ = init_state(M, Fₘ, stateₙ)
function init_state(M::StaticUnivariateZeroState, F, state::StaticUnivariateZeroState)
    init_state(M, F, state.xn0, state.xn1, state.fxn0, state.fxn1)
end


# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::SSecant, F::Callable_Function, x)
    x₀,x₁ = x₀x₁(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# initialize from xs, fxs
function init_state(::SSecant, F, x₀, x₁, fx₀, fx₁)
    SUnivariateZeroState(x₁, x₀, fx₁, fx₀)
end

# return state, convergence_flag
function update_state(method::SSecant, fs, o::SUnivariateZeroState{T,S}, options) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    delta = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)


    if isinf(delta) || isnan(delta)
        return o, true
        #@set! o.message = "Increment `Δx` has issues. "
     end

    @set! o.xn0 = xn1
    @set! o.xn1 -= delta
    @set! o.fxn0 = fxn1
    @set! o.fxn1 = (tmp::S = fs(o.xn1))
    #incfn(o)

    return o, false

end

### has extra Δ
# we store x0,x1,fx0,fx1 and Δ = fx1/f'(x1)
struct SNewtonState{T,S} <: StaticUnivariateZeroState where {T,S}
    xn1::T
    xn0::T
    Δ::T
    fxn1::S
    fxn0::S
end

function init_state(M::SNewton, F::Callable_Function, x)
    x₀ = float(first(x))
    fx₀, Δ = F(x₀)
    x₁ = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end


function init_state(::SNewton, F, x₀, x₁, fx₀, fx₁)
    fx₁, Δ = F(x₁)
    SNewtonState(x₁, x₀, Δ, fx₁, fx₀)
end


function update_state(method::SNewton, fs, o::SNewtonState{T,S}, options) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ = o.Δ


    if isissue(Δ)
        return o, true
    end

    xn0, xn1 = xn1, xn1-Δ
    fxn0 = fxn1
    fxn1, Δ = fs(xn1)
    #incfn(o,2)

    @set! o.xn0 = xn0
    @set! o.xn1 = xn1
    @set! o.Δ = Δ
    @set! o.fxn0 = fxn0
    @set! o.fxn1 = fxn1

    return o, false



end

## Bisection

function init_state(M::SBisection, F::Callable_Function, x)
    x₀, x₁ = float.(_extrema(x))
    fx₀, fx₁ = F(x₀), F(x₁)
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end


function init_state(::SBisection, F, x₀, x₁, fx₀, fx₁)
    assert_bracket(fx₀, fx₁)
    xₘ = Roots._middle(x₀, x₁) # for possibly mixed sign x1, x2
    fxₘ = F(xₘ)
    if sign(fxₘ) * fx₀ < 0
        SUnivariateZeroState(xₘ, x₀, fxₘ, fx₀)
    else
        SUnivariateZeroState(x₁, xₘ, fx₁, fxₘ)
    end
end


function update_state(method::SBisection, F, o::SUnivariateZeroState{T,S}, options) where {T, S}

    x₀, x₁ = o.xn0, o.xn1
    fx₀, fx₁ = o.fxn0, o.fxn1

    xₘ = Roots.__middle(x₀,x₁)
    fxₘ = F(xₘ)

    if sign(fxₘ) * sign(fx₀) < 0
        x₁, fx₁ = xₘ, fxₘ
    else
        x₀, fx₀ = xₘ, fxₘ
    end


    @set! o.xn0 = x₀
    @set! o.xn1 = x₁
    @set! o.fxn0 = fx₀
    @set! o.fxn1 = fx₁

    return o, false



end


## -----

# init_options(M,
#              state::SUnivariateZeroState{T,S},
#              kwargs...
#              ) where {T, S} = init_options(M, T, S; kwargs...)

#     d = kwargs

#     defs = default_tolerances(M, T, S)
#     options = UnivariateZeroOptions(get(d, :xatol, get(d, :xabstol, defs[1])),
#                                     get(d, :xrtol, get(d, :xreltol, defs[2])),
#                                     get(d, :atol,  get(d, :abstol,  defs[3])),
#                                     get(d, :rtol,  get(d, :reltol,  defs[4])),
#                                     get(d, :maxevals,   get(d, :maxsteps, defs[5])),
#                                     get(d, :maxfnevals, defs[6]),
#                                     get(d, :strict,     defs[7]))
#     options
# end

# want abstract type to carry T,S
TS(::SUnivariateZeroState{T,S}) where {T,S} = (T,S)
TS(::SNewtonState{T,S}) where {T,S} = (T,S)
function ffind_zero(f, x, M; kwargs...)

    F = Callable_Function(M, f) # p
    state = init_state(M, F, x)
    T,S = TS(state)
    options = UnivariateZeroOptions(eps(T),eps(T), eps(S), eps(S), 100, 100, false)
    #options = init_options(M, state)
    ffind_zero(M, F, state, options)
end



function ffind_zero(M,F,state,options,l=NullTracks())

#    log_step(l, M, state, :init)
    val, stopped = :not_converged, false
    ctr = 1
    while !stopped
        val, stopped = assess_convergence(M, state, options)
        stopped && break
        ctr > options.maxevals && break
        ctr * fn_argout(M) > options.maxfnevals && break
        state, stopped = update_state(M, F, state, options)
        log_step(l, M, state)
        ctr += 1
    end

    # log_steps
    # log convergence flag (val)

    out = decide_convergence(M,F,state,options, val)



    return out

end

# return convergence flag, stopped
# Convergence flags
# :not_converged
# :x_converged
# :exact_zero
# :f_converged
# :nan
# :inf
function assess_convergence(::Any, state::StaticUnivariateZeroState, options)

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1

    isnan(xn1)

    iszero(fxn1) && return (:exact_zero, true)

    if _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol)
        return (:f_converged, true)
    end

    Δ = abs(xn1 - xn0)
    δ = max(options.xabstol, max(abs(xn1), abs(xn0)) * options.xreltol)
    if Δ ≤ δ
        return (:x_converged, true)
    end

    if isnan(xn1) || isnan(fxn1)
        return (:nan, true)
    end

    if isinf(xn1) || isinf(fxn1)
        return (:inf, true)
    end

    return (:not_converged, false)

end

function decide_convergence(::Any, F, state::StaticUnivariateZeroState, options, val)

    xn0, xn1 = state.xn0, state.xn1
    fxn1 = state.fxn1
    val ∈ (:f_converged, :exact_zero, :converged) && return xn1

    if options.strict
        val == :x_converged && return xn1
    else
        if val ∈ (:x_converged, :not_converged)
            # |fxn| <= relaxed
            _is_f_approx_0(fxn1, xn1, options.abstol, options.reltol, true)
            return xn1
        end
    end

    return NaN*xn1

end
