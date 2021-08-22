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
    find_zero(fs, x0, M, N; p=p, tracks=tracks,verbose=verbose, kwargs...)
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
    state = init_state(M, F, x0)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)

    xstar = find_zero(M, N, F, state, options, l, verbose)
    isnan(xstar) && throw(ConvergenceFailed("Stopped"))

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
                   l::AbstractTracks=NullTracks(),
                   verbose=false
                   )

    incfn(l,2)
    log_step(l, M, state, :init)
    log_method(l,M)
    log_nmethod(l,N)

    quad_ctr = 0
    flag = :not_converged
    ctr = 0
    α = NaN*state.xn1
    while true
        ctr += 1
        flag, converged = assess_convergence(M, state, options)

        converged && break
        ctr >= options.maxevals && break

        state0 = state
        state0, stopped = update_state(M, F, state0, options) # state0 is proposed step

        ## did we find a zero or a bracketing interval?
        if iszero(state0.fxn1)
            state = state0
            break
        elseif sign(state0.fxn0) * sign(state0.fxn1) < 0

            !isa(l, NullTracks) && log_message(l, "Used bracketing method $N on  [$(state0.xn0),$(state0.xn1)], those steps not recorded")

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, state0, Fₙ) # save function calls by using state0 values
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N, Fₙ, stateₙ, optionsₙ))
            break

        end


        ## did we move too far?
        adj = false
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        ts,TB = 1e-3, 1e2 # too small, too big
        if  Δr >= TB * Δx
            adj = true
            r = b + sign(r-b) * TB * Δx  ## too big
        elseif Δr  <= ts *  Δx
            adj = true
            r = b + sign(r - b) * ts * Δx
        end

        @set! state0.xn1 = r
        @set! state0.fxn1 = first(F(r))
        incfn(l)

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0
            a, b = state.xn1, state0.xn1
            fa, fb = state.fxn1, state0.fxn1

            !isa(l, NullTracks) && log_message(l, "Used bracketing method $N on  [$a,$b], those steps not recorded")

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, Fₙ, a, b, fa, fb)
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N,Fₙ,stateₙ,optionsₙ))
            break
        end


        ## did we improve?
        if adj || abs(state0.fxn1) < abs(state.fxn1)
            if isnan(state0.xn1) || isnan(state0.fxn1) || isinf(state0.xn1) || isinf(state0.fxn1)
                break
            end
            state = state0
            log_step(l, M, state)
            quad_ctr = 0
            continue
        end

        ## try quad_vertex, unless that has gotten old
        if quad_ctr > 4
            state = state0
            break
        else
            quad_ctr += 1
            r = quad_vertex(state0.xn1, state0.fxn1, state.xn1, state.fxn1, state.xn0, state.fxn0)

            if isnan(r) || isinf(r)
                state = state0
            else
                fr = F(r)
                incfn(l)

                @set! state0.xn1 = r
                @set! state0.fxn1 = fr
                state = state0
            end
        end
        log_step(l, M, state)
    end

    log_state(l, state)
    verbose && display(l)

    flag, converged = assess_convergence(M, state, options)
    isnan(α) ? decide_convergence(M, F, state, options, flag) : α

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
