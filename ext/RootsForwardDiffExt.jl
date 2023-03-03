module RootsForwardDiffExt

using Roots
import Roots: solve
using ForwardDiff
import ForwardDiff: Dual, derivative, gradient, partials


function solve(ZP::ZeroProblem, M::Roots.AbstractUnivariateZeroMethod, p::Dual{T}; kwargs...) where {T}

    f = ZP.F
    pᵥ = ForwardDiff.value(p)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)

    fₓ = ForwardDiff.partials(f(ForwardDiff.Dual{T}(xᵅ, one(xᵅ)), pᵥ), 1)
    fₚ = ForwardDiff.partials(f(xᵅ, p))

    Dual{T}(xᵅ, - fₚ / fₓ)
end

function solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               p::AbstractArray{<:Dual{T,<:Real}}; kwargs...) where {T}

    f = ZP.F
    pᵥ = map(pᵢ -> pᵢ.value, p)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)

    fₓ = derivative(x -> f(x, pᵥ), xᵅ)
    fₚ = gradient(p -> f(xᵅ, p), pᵥ)
    ∂ = - fₚ / fₓ
    Δs = ntuple(i -> ∂[i] * getindex(partials.(p), i)[i], N)

    Dual{Z}(xᵅ, Δs)
end

end
