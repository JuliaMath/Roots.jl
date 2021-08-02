"""
    Order0()


The `Order0` method is engineered to be a more robust, though possibly
slower, alternative to the other derivative-free root-finding
methods. The implementation roughly follows the algorithm described in
*Personal Calculator Has Key to Solve Any Equation f(x) = 0*, the
SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
The basic idea is to use a secant step. If along the way a bracket is
found, switch to bisection, using `AlefeldPotraShi`.  If the secant
step fails to decrease the function value, a quadratic step is used up
to 4 times.

This is not really 0-order: the secant method has order
1.6...[Wikipedia](https://en.wikipedia.org/wiki/Secant_method#Comparison_with_other_root-finding_methods_
and the the bracketing method has order
1.6180...[Wikipedia](http://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192965-2/S0025-5718-1993-1192965-2.pdf)
so for reasonable starting points and functions, this algorithm should be
superlinear, and relatively robust to non-reasonable starting points.

"""
struct Order0 <: AbstractSecant end

function find_zero(fs, x0, method::Order0;
                   p=nothing,
                   tracks::AbstractTracks=NullTracks(),
                   verbose=false,
                   kwargs...)
    M = Order1()
    N = AlefeldPotraShi()
    F = Callable_Function(M, fs, p)
    find_zero(F, x0, M, N; tracks=tracks,verbose=verbose, kwargs...)
end

function Base.copy(state::UnivariateZeroState{T,S}) where {T, S}
    UnivariateZeroState(state.xn1, state.xn0, state.xstar, copy(state.m),
                        state.fxn1, state.fxn0, state.fxstar, copy(state.fm),
                        state.steps, state.fnevals,
                        state.stopped, state.x_converged,
                        state.f_converged, state.convergence_failed,
                        state.message)
end

function Base.copy!(dst::UnivariateZeroState{T,S}, src::UnivariateZeroState{T,S}) where {T,S}
    dst.xn1 = src.xn1
    dst.xn0 = src.xn0
    empty!(dst.m); append!(dst.m, src.m)
    dst.fxn1 = src.fxn1
    dst.fxn0 = src.fxn0
    empty!(dst.fm); append!(dst.fm, src.fm)
    dst.steps = src.steps
    dst.fnevals = src.fnevals
    dst.stopped = src.stopped
    dst.x_converged = src.x_converged
    dst.f_converged = src.f_converged
    dst.convergence_failed = src.convergence_failed
    dst.message = src.message
    nothing
end

# todo: consolidate this with find_zero(M,N,f, x0)...
function find_zero(fs, x0,
                   M::AbstractUnivariateZeroMethod,
                   N::AbstractBracketing;
                   p = nothing,
                   tracks::AbstractTracks=NullTracks(),
                   verbose=false,
                   kwargs...)

    F = Callable_Function(M, fs, p)
    #    state = init_state(M, F, x0)
    state = init_state(M, F, x0)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)

    xstar = find_zero(M, N, F, state, options, l)
    verbose &&  show_trace(M, N, state, l)

    isnan(xstar) && throw(ConvergenceFailed("Stopped at: xn = $(state.xn1). $(state.message)"))

    return xstar

end


# Robust version using some tricks: idea from algorithm described in
# [The SOLVE button from the
# HP-34]C(http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
# * use bracketing method if one identifed
# * limit steps so as not too far or too near the previous one
# * if not decreasing, use a quad step upto 4 times to bounce out of trap, if possible
# First uses M, then N if bracket is identified
function find_zero(M::AbstractUnivariateZeroMethod,
                   N::AbstractBracketing,
                   F,
                   state::AbstractUnivariateZeroState,
                   options::UnivariateZeroOptions,
                   l::AbstractTracks=NullTracks()
                   )


    log_step(l, M, state, :init)
    state0 = copy(state)
    quad_ctr = 0

    while true

        if assess_convergence(M, state, options)
            break
        end

        copy!(state0, state)
        update_state(M, F, state0, options) # state0 is proposed step

        adj = false
        ## did we find a zero or a bracketing interval?
        if iszero(state0.fxn1)
            copy!(state, state0)
            state.xstar, state.fxstar = state.xn1, state.fxn1
            state.f_converged = true
            break
        elseif sign(state0.fxn0) * sign(state0.fxn1) < 0
            copy!(state, state0)
            a, b = state0.xn0, state0.xn1 # could save some fn calls here
            run_bisection(N, F, (a, b), state)
            break
        end


        ## did we move too far?
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        ts,TB = 1e-3, 1e2 # too small, too big
        if  Δr >= TB * Δx
            adj = true
            r = b + sign(r-b) * TB * Δx  ## too big
            state0.xn1 = r
            state0.fxn1 = F(r)
            incfn(state)
        elseif Δr  <= ts *  Δx
            adj = true
            r = b + sign(r - b) * ts * Δx
            fr = first(F(r))
            incfn(state)
            state0.xn1 = r
            state0.fxn1 = fr
        end

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0
            state.xn0, state.fxn0 = state.xn1, state.fxn1
            state.xn1, state.fxn1 = state0.xn1, state0.fxn1
            a, b = state.xn0, state.xn1
            run_bisection(N, F, (a, b), state)
            break
        end


        ## did we improve?
        if adj || abs(state0.fxn1) < abs(state.fxn1)
            if isnan(state0.xn1) || isnan(state0.fxn1) || isinf(state0.xn1) || isinf(state0.fxn1)
                break
            end
            copy!(state, state0)
            log_step(l, M, state)
            incsteps(state)
            quad_ctr = 0
            continue
        end

        ## try quad_vertex, unless that has gotten old
        if quad_ctr > 4
            copy!(state, state0)
            state.stopped = true
            break
        else
            quad_ctr += 1
            r = quad_vertex(state0.xn1, state0.fxn1, state.xn1, state.fxn1, state.xn0, state.fxn0)

            if isnan(r) || isinf(r)
                copy!(state, state0)
            else

                fr = F(r)
                incfn(state)

                state0.xn1 = r
                state0.fxn1 = fr
                incfn(state)
                copy!(state, state0)
            end
        end

        log_step(l, M, state)
        incsteps(state)
    end

    decide_convergence(M, F, state, options)
end

# Switch to bracketing method
function run_bisection(N::AbstractBracketing, f, ab, state)
    steps, fnevals = state.steps, state.fnevals
    f = Callable_Function(N, f)
    init_state!(state, N, f; clear=true)
    find_zero(N, f, state, init_options(N, state))
    a, b = _extrema(ab)
    u,v = a > b ? (b, a) : (a, b)
    state.message *= "Bracketing used over ($u, $v), those steps not shown. "
    return nothing
end
