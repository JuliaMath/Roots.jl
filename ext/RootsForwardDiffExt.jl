module RootsForwardDiffExt

using Roots
import Roots: solve
using ForwardDiff
import ForwardDiff: Dual, value, partials

function ∂f∂p(ZP, M, p, T; kwargs...)
    f = ZP.F
    pᵥ = map(value, p)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)
    𝐱ᵅ = Dual{T}(xᵅ, one(xᵅ))

    fₓ = partials(f(𝐱ᵅ, pᵥ), 1)  # derivative(x -> f(x, pᵥ), xᵅ)
    fₚ = partials(f(xᵅ, p))

    Dual{T}(xᵅ, - fₚ / fₓ)
end


function solve(ZP::ZeroProblem, M::Roots.AbstractUnivariateZeroMethod, p::Dual{T}; kwargs...) where {T}

    ∂f∂p(ZP, M, p, T; kwargs...)

end

function solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               p::AbstractArray{<:Dual{T,<:Real}}; kwargs...) where {T}
    ∂f∂p(ZP, M, p, T; kwargs...)
end

end
