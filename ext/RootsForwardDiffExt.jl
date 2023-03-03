module RootsForwardDiffExt

using Roots
import Roots: solve
using ForwardDiff
import ForwardDiff: Dual, value, partials

function ∂f∂p(ZP, M, 𝐩, T; kwargs...)
    f = ZP.F
    pᵥ = value.(𝐩)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)
    𝐱ᵅ = Dual{T}(xᵅ, one(xᵅ))

    fₓ = partials(f(𝐱ᵅ, pᵥ), 1)
    fₚ = partials(f(xᵅ, 𝐩))

    Dual{T}(xᵅ, - fₚ / fₓ)
end

function solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               𝐩::Dual{T}; kwargs...) where {T}
    ∂f∂p(ZP, M, 𝐩, T; kwargs...)
end

function solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               𝐩::AbstractArray{<:Dual{T,<:Real}}; kwargs...) where {T}
    ∂f∂p(ZP, M, 𝐩, T; kwargs...)
end

end
