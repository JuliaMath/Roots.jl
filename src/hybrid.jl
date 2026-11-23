## Init for hybrid method -- start with a non-bracketing, finish with bracketing

## When passing 2 methods, any parameters must be passed as a named argument through
## the keyword p
function init(
    𝑭𝑿::ZeroProblem,
    M::AbstractNonBracketingMethod,
    N::AbstractBracketingMethod,
    p′=nothing;
    p=nothing,
    tracks=NullTracks(),
    kwargs...,
)
    F = Callable_Function(M, 𝑭𝑿.F, something(p′, p, missing))
    state = init_state(M, F, 𝑭𝑿.x₀)
    options = init_options(M, state; kwargs...)
    incfn(tracks, initial_fncalls(M))
    ZeroProblemIterator(M, N, F, state, options, tracks)
end

# Robust version using some tricks: idea from algorithm described in
# [The SOLVE button from the
# HP-34C](https://hparchive.com/Journals/HPJ-1979-12.pdf).
# * use bracketing method if one is identified
# * limit steps to at most 100 times the previous one (nearly horizontal secants)
# * if a step rounds to zero away from a root, take a small step instead of stopping
# * if not decreasing, use a quad step upto 4 times to bounce out of trap, if possible
# First uses M, then N if bracket is identified
function solve!(𝐙::ZeroProblemIterator{𝐌,𝐍};) where {𝐌,𝐍<:AbstractBracketingMethod}
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
        crossed = sign(state0.fxn0) * sign(state0.fxn1) < 0
        if iszero(state0.fxn1) || (crossed && is_approx_zero_f(M, state0, options))
            # an exact zero, or crossed the zero within the tolerance of M; a
            # non-crossing step is only checked after limiting its size, as a
            # far-off step may pass the relative tolerance with |f| not small
            state = state0
            break
        elseif crossed
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

            return α
        end

        ## Kahan's two amendments to the step
        adj = false
        r, a, b = state0.xn1, state.xn0, state.xn1
        Δr = abs(r - b)
        Δx = abs(b - a)
        if Δr >= 100 * Δx
            # too big: limit nearly horizontal secants (Kahan: |γ - β| ≤ 100|β - α|)
            adj = true
            r = b + sign(r - b) * 100 * Δx
        elseif iszero(Δr) && !iszero(Δx) && !is_approx_zero_f(M, state, options, true)
            # step rounded to zero away from a root: persevere (Kahan, Fig. 7)
            adj = true
            r = b + sign(b - a) * Δx / 1000
        end

        if adj
            @reset state0.xn1 = r
            @reset state0.fxn1 = first(F(r))
            incfn(l)
        end

        # a sign change after shortening?
        if sign(state.fxn1) * sign(state0.fxn1) < 0 && !is_approx_zero_f(M, state0, options)
            log_step(l, M, state)
            a, b = state.xn1, state0.xn1
            fa, fb = state.fxn1, state0.fxn1
            !isa(l, NullTracks) && log_message(l, "Used bracketing method $N on  [$a,$b]")

            Fₙ = Callable_Function(N, F)
            stateₙ = init_state(N, Fₙ, a, b, fa, fb)
            optionsₙ = init_options(N, stateₙ)
            α = solve!(init(N, Fₙ, stateₙ, optionsₙ, l))

            log_method(l, M)

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
                fr = first(F(r))
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

    isnan(α) ? decide_convergence(M, F, state, options, flag) : α
end

function find_zero(
    fs,
    x0,
    M::AbstractUnivariateZeroMethod,
    N::AbstractBracketingMethod,
    p′=nothing;
    kwargs...,
)
    𝐏 = ZeroProblem(fs, x0)
    solve!(init(𝐏, M, N, p′; kwargs...))
end
