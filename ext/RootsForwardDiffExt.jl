module RootsForwardDiffExt

using Roots
import Roots: solve
using ForwardDiff
import ForwardDiff: Dual, derivative, gradient, partials


function solve(ZP::ZeroProblem, M::Roots.AbstractUnivariateZeroMethod, p::Dual{T}; kwargs...) where {T}

    f = ZP.F
    pᵥ = ForwardDiff.value(p)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)

    fₓ = derivative(x -> f(x, pᵥ), xᵅ)
    fₚ = derivative(p -> f(xᵅ, p), pᵥ)

    Dual{T}(xᵅ, - fₚ / fₓ)
end

function solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               p::AbstractVector{<:Dual{Z,T,N}}; kwargs...) where {Z,T<:Real,N}

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
