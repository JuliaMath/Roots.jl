## Init for hybrid method -- start with a non-bracketing, finish with bracketing

## When passing 2 methods, any parameters must be passed as a named argument through
## the keyword p
function init(
    𝑭𝑿::ZeroProblem,
    M::AbstractNonBracketingMethod,
    N::AbstractBracketingMethod,
    p′=nothing;
    p=nothing,
    verbose::Bool=false,
    tracks=NullTracks(),
    kwargs...,
)
    F = Callable_Function(M, 𝑭𝑿.F, something(p′, p, missing))
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    l =  (verbose && isa(tracks, NullTracks)) ? Tracks() : tracks
    incfn(l, initial_fncalls(M))
    ZeroProblemIterator(M, N, F, state, options, l)
end

# Robust version using some tricks: idea from algorithm described in
# [The SOLVE button from the
# HP-34]C(http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).
# * use bracketing method if one is identified
# * limit steps so as not too far or too near the previous one
# * if not decreasing, use a quad step upto 4 times to bounce out of trap, if possible
# First uses M, then N if bracket is identified
function solve!(
    𝐙::ZeroProblemIterator{𝐌,𝐍};
    verbose=false,
) where {𝐌,𝐍<:AbstractBracketingMethod}
    M, N, F, state, options, l = 𝐙.M, 𝐙.N, 𝐙.F, 𝐙.state, 𝐙.options, 𝐙.logger

    incfn(l, 2)
    log_step(l, M, state; init=true)

    log_method(l, M)
    log_nmethod(l, N)

    quad_ctr = 0
    flag = :not_converged
    ctr = 0
    α = nan(typeof(state.xn1)) * state.xn1
    while true
        ctr += 1
        flag, converged = assess_convergence(M, state, options)

        converged && break
        ctr >= options.maxiters && break

        state0 = state
        state0, stopped = update_state(M, F, state0, options) # state0 is proposed step

        ## did we find a zero or a bracketing interval?
        if iszero(state0.fxn1)
            state = stagte0
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

            log_method(l, M)
            verbose && display(l)

            return α
        end

        ## did we move too far?
        adj = false
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        ts, TB = one(r) / 1000, 100 * one(r) # too small, too big
        if Δr >= TB * Δx
            adj = true
            r = b + sign(r - b) * TB * Δx  ## too big
        elseif Δr <= ts * Δx
            adj = true
            r = b + sign(r - b) * ts * Δx
        end

        @reset state0.xn1 = r
        @reset state0.fxn1 = first(F(r))
        incfn(l)

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0
            log_step(l, M, state)
            a, b = state.xn1, state0.xn1
            fa, fb = state.fxn1, state0.fxn1
            !isa(l, NullTracks) && log_message(l, "Used bracketing method $N on  [$a,$b]")

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, Fₙ, a, b, fa, fb)
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N, Fₙ, stateₙ, optionsₙ, l))

            log_method(l, M)
            verbose && display(l)

            return α
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

                @reset state0.xn1 = r
                @reset state0.fxn1 = fr
                state = state0
            end
        end
        log_step(l, M, state)
    end

    val, stopped = assess_convergence(M, state, options)
    α = decide_convergence(M, F, state, options, val)

    log_convergence(l, val)
    log_last(l, α)
    verbose && display(l)

    isnan(α) ? decide_convergence(M, F, state, options, flag) : α
end

function find_zero(
    fs,
    x0,
    M::AbstractUnivariateZeroMethod,
    N::AbstractBracketingMethod,
    p′=nothing;
    verbose=false,
    kwargs...,
)
    𝐏 = ZeroProblem(fs, x0)
    solve!(init(𝐏, M, N, p′; verbose=verbose, kwargs...), verbose=verbose)
end
