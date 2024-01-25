# View find_zero as solving `f(x, p) = 0` for `xᵅ(p)`.
# This is implicitly defined. By the implicit function theorem, we have:
# ∇f = 0 => ∂/∂ₓ f(xᵅ, p) ⋅ ∂xᵅ/∂ₚ + ∂/∂ₚf(x\^α, p) ⋅ I = 0
# or ∂xᵅ/∂ₚ = - ∂/∂ₚ f(xᵅ, p)  / ∂/∂ₓ f(xᵅ, p)

function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, _, _, Δp),
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    p;
    kwargs...,
)
    @show :frule, p
    xᵅ = solve(ZP, M, p; kwargs...)

    # Use a single reverse-mode AD call with `rrule_via_ad` if `config` supports it?
    F = p -> Callable_Function(M, ZP.F, p)
    fₓ(x) = first(F(p)(x))
    fₚ(p) = first(F(p)(xᵅ))
    fx = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), true), fₓ, xᵅ)[2]
    fp = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), Δp), fₚ, p)[2]

    xᵅ, -fp / fx
end



function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, _, _, Δp),
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    ::Nothing;
    kwargs...,
)

    @show :frule, :nothing
    xᵅ = solve(ZP, M, p; kwargs...)

    # Use a single reverse-mode AD call with `rrule_via_ad` if `config` supports it?
    F = p -> Callable_Function(M, ZP.F, p)
    fₓ(x) = first(F(p)(x))
    fₚ(p) = first(F(p)(xᵅ))
    fx = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), true), fₓ, xᵅ)[2]
    fp = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), Δp), fₚ, p)[2]

    xᵅ, -fp / fx
end

## modified from
## https://github.com/gdalle/ImplicitDifferentiation.jl/blob/main/src/implicit_function.jl
function ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    p;
    kwargs...,
)
    @show :rrule, p
    xᵅ = solve(ZP, M, p; kwargs...)

    f(x, p) = first(Callable_Function(M, ZP.F, p)(x))
    _, pullback_f = ChainRulesCore.rrule_via_ad(rc, f, xᵅ, p)
    _, fx, fp = pullback_f(true)
    yp = -fp / fx

    function pullback_solve_ZeroProblem(dy)
        dp = yp * dy
        return (
            ChainRulesCore.NoTangent(),
            ChainRulesCore.NoTangent(),
            ChainRulesCore.NoTangent(),
            dp,
        )
    end

    return xᵅ, pullback_solve_ZeroProblem
end

function ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    ::Nothing;
    kwargs...,
)
    @show :rrule, nothing
    xᵅ = solve(ZP, M; kwargs...)
    𝐹 = typeof(ZP.F)
    @show 𝐹
    f(x, p) = first(Callable_Function(M, p)(x))
    _, pullback_f = ChainRulesCore.rrule_via_ad(rc, f, xᵅ, ZP.F)
    _, fx, fp = pullback_f(true)
    @show fx, first(fp)
    yp = -first(fp) / fx

    function pullback_solve_ZeroProblem(dy)
        dp = yp * dy
        @show dp
        return (
            ChainRulesCore.NoTangent(),
            ChainRulesCore.NoTangent(),
            fx,#ChainRulesCore.NoTangent(),
            dp
        )
    end

    return xᵅ, pullback_solve_ZeroProblem
end
