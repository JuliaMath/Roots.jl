module RootsForwardDiffExt

using Roots
using ForwardDiff
import ForwardDiff: Dual, value, partials, Partials, derivative, gradient!

# What works
# F(p) = find_zero(f, x0, M, p)
# G(p) = find_zero(𝐺(p), x0, M)
#                           F                    G
# ForwardDiff.derivative    ✓                    x (wrong answer, 0.0)
# ForwardDiff.gradient      ✓                    x (wrong answer, 0.0)
# ForwardDiff.hessian       ✓                    x (wrong answer, 0.0)
# Zygote.gradient           ✓                    ✓
# Zygote.hessian            ✓                    x (wrong answer!)
# Zygote.hessian_reverse    ✓                    x (MethodError)

function Roots.solve(ZP::ZeroProblem,
                     M::Roots.AbstractUnivariateZeroMethod,
                     𝐩::Dual{T};
                     kwargs...) where {T}


    # p_and_dp = 𝐩
    p, dp = value.(𝐩), partials.(𝐩)

    xᵅ = solve(ZP, M, p; kwargs...)

    f(x, p) = first(Roots.Callable_Function(M, ZP.F, p)(x))
    fₓ = derivative(_x -> f(_x, p), xᵅ)
    fₚ = derivative(_p -> f(xᵅ, _p), p)

    # x and dx
    dx = - (fₚ * dp) / fₓ

    Dual{T}(xᵅ, dx)
end

# cf https://discourse.julialang.org/t/custom-rule-for-differentiating-through-newton-solver-using-forwarddiff-works-for-gradient-fails-for-hessian/93002/22
function Roots.solve(ZP::ZeroProblem,
                     M::Roots.AbstractUnivariateZeroMethod,
                     𝐩::AbstractArray{<:Dual{T,R,N}};
                     kwargs...) where {T,R,N}


    # p_and_dp = 𝐩
    p, dp = value.(𝐩), partials.(𝐩)
    xᵅ = solve(ZP, M, p; kwargs...)

    f(x, p) = first(Roots.Callable_Function(M, ZP.F, p)(x))
    fₓ = derivative(_x -> f(_x, p), xᵅ)
    fₚ = similar(p)  # <-- need this, not output of gradient(p->f(x,p), p)
    gradient!(fₚ, _p -> f(xᵅ, _p), p)

    # x_and_dx
    dx = - (fₚ' * dp) / fₓ

    Dual{T}(xᵅ, Partials(ntuple(k -> dx[k], Val(N))))
end
end
