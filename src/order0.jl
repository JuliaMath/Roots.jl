"""
    Order0()


The `Order0` method is engineered to be a more robust, though possibly
slower, alternative to the other derivative-free root-finding
methods. The implementation roughly follows the algorithm described in
*Personal Calculator Has Key to Solve Any Equation ``f(x) = 0``*, the
SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
The basic idea is to use a secant step. If along the way a bracket is
found, switch to a bracketing algorithm, using `AlefeldPotraShi`.  If the secant
step fails to decrease the function value, a quadratic step is used up
to ``4`` times.

This is not really ``0``-order: the secant method has order
``1.6...`` [Wikipedia](https://en.wikipedia.org/wiki/Secant_method#Comparison_with_other_root-finding_methods)
and the the bracketing method has order
``1.6180...`` [Wikipedia](http://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192965-2/S0025-5718-1993-1192965-2.pdf)
so for reasonable starting points and functions, this algorithm should be
superlinear, and relatively robust to non-reasonable starting points.

"""
struct Order0 <: AbstractSecantMethod end

# special case Order0 to be hybrid
function init(
    𝑭𝑿::ZeroProblem,
    M::Order0,
    p′=nothing;
    p=nothing,
    verbose::Bool=false,
    tracks=NullTracks(),
    kwargs...,
)
    p = p′ === nothing ? p : p′
    init(𝑭𝑿, Secant(), AlefeldPotraShi(); p=p, verbose=verbose, tracks=tracks, kwargs...)
end

## When passing 2 methods, any parameters must be passed as a named argument through
## the keyword p
function init(
    𝑭𝑿::ZeroProblem,
    M::AbstractNonBracketingMethod,
    N::AbstractBracketingMethod;
    p=nothing,
    verbose::Bool=false,
    tracks=NullTracks(),
    kwargs...,
)
    F = Callable_Function(M, 𝑭𝑿.F, p)
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    l = Tracks(verbose, tracks, state)
    incfn(l, initial_fncalls(M))
    ZeroProblemIterator(M, N, F, state, options, l)
end

# Robust version using some tricks: idea from algorithm described in
# [The SOLVE button from the
# HP-34]C(http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
# * use bracketing method if one identifed
# * limit steps so as not too far or too near the previous one
# * if not decreasing, use a quad step upto 4 times to bounce out of trap, if possible
# First uses M, then N if bracket is identified
function solve!(𝐙::ZeroProblemIterator{𝐌,𝐍}; verbose=false) where {𝐌,𝐍<:AbstractBracketingMethod}
    M, N, F, state, options, l = 𝐙.M, 𝐙.N, 𝐙.F, 𝐙.state, 𝐙.options, 𝐙.logger

    incfn(l, 2)
    log_step(l, M, state; init=true)

    quad_ctr = 0
    flag = :not_converged
    ctr = 0
    α = NaN * state.xn1
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
            log_step(l, M, state0)
            !isa(l, NullTracks) && log_message(
                l,
                "Used bracketing method $N on  [$(min(state0.xn0,state0.xn1)),$(max(state0.xn0,state0.xn1))]",
            )

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, state0, Fₙ) # save function calls by using state0 values
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N, Fₙ, stateₙ, optionsₙ, l))
            break
        end

        ## did we move too far?
        adj = false
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        ts, TB = 1e-3, 1e2 # too small, too big
        if Δr >= TB * Δx
            adj = true
            r = b + sign(r - b) * TB * Δx  ## too big
        elseif Δr <= ts * Δx
            adj = true
            r = b + sign(r - b) * ts * Δx
        end

        @set! state0.xn1 = r
        @set! state0.fxn1 = first(F(r))
        incfn(l)

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0
            log_step(l, M, state)
            a, b = state.xn1, state0.xn1
            fa, fb = state.fxn1, state0.fxn1
            !isa(l, NullTracks) && log_message(
                l,
                "Used bracketing method $N on  [$a,$b]",
            )

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, Fₙ, a, b, fa, fb)
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N, Fₙ, stateₙ, optionsₙ, l))
            break
        end

        ## did we improve?
        if adj || abs(state0.fxn1) < abs(state.fxn1)
            if isnan(state0.xn1) ||
               isnan(state0.fxn1) ||
               isinf(state0.xn1) ||
               isinf(state0.fxn1)
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
            r = quad_vertex(
                state0.xn1,
                state0.fxn1,
                state.xn1,
                state.fxn1,
                state.xn0,
                state.fxn0,
            )

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


    val, stopped = assess_convergence(M, state, options)
    log_convergence(l, val)
    log_method(l, M)
    log_nmethod(l, N)
    log_last(l, α)
    verbose && display(l)

    isnan(α) ? decide_convergence(M, F, state, options, flag) : α
end

function find_zero(
    fs,
    x0,
    M::AbstractUnivariateZeroMethod,
    N::AbstractBracketingMethod;
    verbose=false,
    kwargs...,
)
    𝐏 = ZeroProblem(fs, x0)
    solve!(init(𝐏, M, N; verbose=verbose, kwargs...), verbose=verbose)
end

# Switch to bracketing method
# deprecate soon, not used
function run_bisection(N::AbstractBracketingMethod, f, ab, state)
    steps, fnevals = state.steps, state.fnevals
    f = Callable_Function(N, f)
    init_state!(state, N, f; clear=true)
    find_zero(N, f, state, init_options(N, state))
    a, b = _extrema(ab)
    u, v = a > b ? (b, a) : (a, b)
    state.message *= "Bracketing used over ($u, $v), those steps not shown. "
    return nothing
end
