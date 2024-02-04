# View find_zero as solving `f(x, p) = 0` for `xᵅ(p)`.
# This is implicitly defined. By the implicit function theorem, we have:
# ∇f = 0 => ∂/∂ₓ f(xᵅ, p) ⋅ ∂xᵅ/∂ₚ + ∂/∂ₚf(x\^α, p) ⋅ I = 0
# or ∂xᵅ/∂ₚ = - ∂/∂ₚ f(xᵅ, p)  / ∂/∂ₓ f(xᵅ, p)

# There are two cases considered
# F(p) = find_zero(f(x,p), x₀, M, p) # f a function
# G(p) = find_zero(𝐺(p), x₀, M)      # 𝐺 a functor
# For G(p) first order derivatives are working
# **but** hessian is not with Zygote. *MOREOVER* it fails
# with the **wrong answer** not an error.
#
# (`Zygote.hessian` calls `ForwardDiff` and that isn't working with a functor;
# `Zygote.hessian_reverse` doesn't seem to work here, though perhaps
# that is fixable.)


# this assumes a function and a parameter `p` passed in
import ChainRulesCore: Tangent, NoTangent, frule, rrule
function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, _, _, Δp),
    ::typeof(solve),
    ZP::ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod,
    p;
    kwargs...,
)
    xᵅ = solve(ZP, M, p; kwargs...)

    # Use a single reverse-mode AD call with `rrule_via_ad` if `config` supports it?
    F = p -> Roots.Callable_Function(M, ZP.F, p)
    fₓ(x) = first(F(p)(x))
    fₚ(p) = first(F(p)(xᵅ))
    fx = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), true), fₓ, xᵅ)[2]
    fp = ChainRulesCore.frule_via_ad(config, (ChainRulesCore.NoTangent(), Δp), fₚ, p)[2]

    xᵅ, -fp / fx
end

# Case of Functor carrying parameters
ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    xdots,
    ::typeof(solve),
    ZP::Roots.ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod,
    ::Nothing;
    kwargs...,
) =
    frule(config, xdots, solve, ZP, M; kwargs...)

function ChainRulesCore.frule(
    config::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasForwardsMode},
    (_, Δq, _),
    ::typeof(solve),
    ZP::Roots.ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod;
    kwargs...,
)
    # no `p`; make ZP.F the parameter (issue 408)
    foo = ZP.F
    zprob2 = ZeroProblem(|>, ZP.x₀)
    nms = fieldnames(typeof(foo))
    nt = NamedTuple{nms}(getfield(foo, n) for n ∈ nms)
    dfoo = Tangent{typeof(foo)}(;nt...)

    return frule(config,
                 (NoTangent(), NoTangent(), NoTangent(), dfoo),
                 Roots.solve, zprob2, M, foo)
end


##

## modified from
## https://github.com/gdalle/ImplicitDifferentiation.jl/blob/main/src/implicit_function.jl
# this is for passing a parameter `p`
function ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
    ::typeof(solve),
    ZP::ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod,
    p;
    kwargs...,
)
    xᵅ = solve(ZP, M, p; kwargs...)

    f(x, p) = first(Roots.Callable_Function(M, ZP.F, p)(x))
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

# this assumes a functor 𝐺(p) for the function *and* no parameter
ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
    ::typeof(solve),
    ZP::ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod,
    ::Nothing;
    kwargs...,
) =
    ChainRulesCore.rrule(rc, solve, ZP, M; kwargs...)


function ChainRulesCore.rrule(
    rc::ChainRulesCore.RuleConfig{>:ChainRulesCore.HasReverseMode},
    ::typeof(solve),
    ZP::ZeroProblem,
    M::Roots.AbstractUnivariateZeroMethod;
    kwargs...,
)


    𝑍𝑃 = ZeroProblem(|>, ZP.x₀)
    xᵅ = solve(ZP, M; kwargs...)
    f(x, p) = first(Roots.Callable_Function(M, 𝑍𝑃.F, p)(x))

    _, pullback_f = ChainRulesCore.rrule_via_ad(rc, f, xᵅ, ZP.F)
    _, fx, fp = pullback_f(true)

    yp = NamedTuple{keys(fp)}(-fₚ/fx for fₚ ∈ values(fp))

    function pullback_solve_ZeroProblem(dy)
        dF = ChainRulesCore.Tangent{typeof(ZP.F)}(; yp...)

        dZP = ChainRulesCore.Tangent{typeof(ZP)}(;
                                                 F = dF,
                                                 x₀ = ChainRulesCore.NoTangent()
                                                 )

        dsolve = ChainRulesCore.NoTangent()
        dM = ChainRulesCore.NoTangent()
        dp = ChainRulesCore.NoTangent()

        return dsolve, dZP, dM, dp
    end

    return xᵅ, pullback_solve_ZeroProblem
end
