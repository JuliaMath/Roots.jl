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
    @show ZP.F
    @show fp
    @show fx
    @show yp
    function pullback_solve_ZeroProblem(dy)
        @show dy
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
    @show :rrule, nothing, :XXX_not_working_for_hessian
    𝑓(x,p) = p(x)
    𝑍𝑃 = ZeroProblem(𝑓, ZP.x₀)
    @show ZP.F

    xᵅ = solve(ZP, M; kwargs...)
    f(x, p) = first(Callable_Function(M, 𝑍𝑃.F, p)(x))
    # XXX --> this is failing
    _, pullback_f = ChainRulesCore.rrule_via_ad(rc, f, xᵅ, ZP.F)
    _, fx, fp = pullback_f(true)
    @show isa(first(fp), Number), typeof(first(fp))
    yp = NamedTuple{keys(fp)}(values(fp)./(-fx))
    @show [getfield(ZP.F,k) for k in keys(fp)]
    @show fp
    @show fx
    @show yp
    function pullback_solve_ZeroProblem(dy)
        #𝐹 = hasproperty(ZP.F, :backing) ? getfield(ZP.F, :backing) : ZP.F
        #dys = [getfield(𝐹,k) / fx for k in keys(fp)]
        #@show dys
        #dys = [getfield(fp,k)*getfield(ZP.F,k) / fx for k in keys(fp)]
        #@show dy
        dF = ChainRulesCore.Tangent{typeof(ZP.F)}(; yp...) # but not right for hessian
        dsolve = ChainRulesCore.NoTangent()

        dZP = ChainRulesCore.Tangent{typeof(ZP)}(;
                                                 F = dF,
                                                 x₀ = ChainRulesCore.NoTangent()
                                                 )

        dM = ChainRulesCore.NoTangent()
        dp = ChainRulesCore.NoTangent()
        dF = ChainRulesCore.@thunk(values(yp) .* dy)
        return dsolve, dZP, dM, dp, dF
    end

    return xᵅ, pullback_solve_ZeroProblem
end
