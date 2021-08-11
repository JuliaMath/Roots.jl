using Setfield


struct SSecant <: AbstractSecant end
export SSecant
fn_argout(::SSecant) = 1

struct SUnivariateZeroState{T,S} <: AbstractUnivariateZeroState where {T,S}
    xn1::T
    xn0::T
    xstar::T
    m::Vector{T}
    fxn1::S
    fxn0::S
    fxstar::S
    fm::Vector{S}
    steps::Int
    fnevals::Int
    stopped::Bool             # stopped, butmay not have converged
    x_converged::Bool         # converged via |x_n - x_{n-1}| < ϵ
    f_converged::Bool         # converged via |f(x_n)| < ϵ
    convergence_failed::Bool
    message::String
end

function _init_state(M::SSecant,
    x₀::T, x₁::T, fx₀::S, fx₁::S;
    m=T[],
    fm=S[],
    xstar = nan(T)*oneunit(T),
    fxstar = nan(S)*oneunit(S),
    steps = 0,
    fnevals= 0,
    stopped=false,
    x_converged=false, f_converged=false,
    convergence_failed=false,
    message::String="") where {T,S}

    UnivariateZeroState(x₁, x₀, xstar, m,
                        fx₁, fx₀, fxstar, fm,
                        steps, fnevals,
                        stopped, x_converged, f_converged,
                        convergence_failed,
                        message)

end


function update_state(method::SSecant, fs, o::UnivariateZeroState{T,S}, options) where {T, S}

    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1

    delta = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)


    if isinf(delta) || isnan(delta)
        o.stopped = true
        @set! o.message = "Increment `Δx` has issues. "
        return o
     end

    @set! o.xn0 = xn1
    @set! o.xn1 -= delta
    @set! o.fxn0 = fxn1
    @set! o.fxn1 = (tmp::S = fs(o.xn1))
    incfn(o)

    return o

end


function ffind_zero(M,F,state,options,l=NullTracks())

    log_step(l, M, state, :init)

    while !assess_convergence(M, state, options)
        state = update_state(M, F, state, options)
        log_step(l, M, state)
        incsteps(state)
    end

    decide_convergence(M,F,state,options)
    return state

xBgXU
end
