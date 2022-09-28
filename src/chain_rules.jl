# View find_zero as solving `f(x, p) = 0` for `xᵅ(p)`.
# This is implicitly defined. By the implicit function theorem, we have:
# ∇f = 0 => ∂/∂ₓ f(xᵅ, p) ⋅ ∂xᵅ/∂ₚ + ∂/∂ₚf(x\^α, p) ⋅ I = 0
# or ∂xᵅ/∂ₚ =  ∂/∂ₚf(xᵅ, p)  / ∂/∂ₓ f(xᵅ, p)

# does this work?
# It doesn't pass a few of the tests of ChainRulesTestUtils
function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig,
    (Δself, Δp),
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    p;
    kwargs...)


    xᵅ = solve(ZP, M, p; kwargs...)

    F = p -> Callable_Function(M, ZP.F, p)
    fₓ(x) = first(F(p)(x))
    fₚ(p) = - first(F(p)(xᵅ))

    function pushforward_find_zero(fₓ, fₚ, xᵅ, p, Δp)
        # is scalar?
        o = typeof(p) == eltype(p) ? one(p) : ones(eltype(p), size(p))
        fx = ChainRulesCore.frule_via_ad(config,
                                         (ChainRulesCore.NoTangent(), o),
                                         fₓ, xᵅ)[2]
        fp = ChainRulesCore.frule_via_ad(config,
                                         (ChainRulesCore.NoTangent(), o),
                                         fₚ, p)[2]

        dp = ChainRulesCore.unthunk(Δp)
        δ = - (fp * dp) / fx
        δ
    end

    xᵅ, pushforward_find_zero(fₓ, fₚ, xᵅ, p, Δp)

end

## modified from
## https://github.com/gdalle/ImplicitDifferentiation.jl/blob/main/src/implicit_function.jl
function ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig,
    ::typeof(solve),
    ZP::ZeroProblem,
    M::AbstractUnivariateZeroMethod,
    p;
    kwargs...)

    xᵅ = solve(ZP, M, p; kwargs...)
    F = p -> Callable_Function(M, ZP.F, p)
    fₓ(x) = first(F(p)(x))
    fₚ(p) = - first(F(p)(xᵅ))

    pullback_Aᵀ = last ∘ ChainRulesCore.rrule_via_ad(rc, fₓ, xᵅ)[2]
    pullback_Bᵀ = last ∘ ChainRulesCore.rrule_via_ad(rc, fₚ, p)[2]

    function pullback_find_zero(dy)
        dy =  ChainRulesCore.unthunk(dy)
        u = inv(pullback_Aᵀ(1/dy))
        dx = pullback_Bᵀ(u)
        return (ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(),
                ChainRulesCore.NoTangent(), dx)
    end

    return xᵅ, pullback_find_zero
end
