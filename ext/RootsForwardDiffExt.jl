module RootsForwardDiffExt

using Roots
using ForwardDiff
import ForwardDiff: Dual, value, partials

# For ForwardDiff we add a `solve` method for Dual types
function Roots.solve(ZP::ZeroProblem,
               M::Roots.AbstractUnivariateZeroMethod,
               𝐩::Union{Dual{T},
                        AbstractArray{<:Dual{T,<:Real}}
                        };
               kwargs...) where {T}
    f = ZP.F
    pᵥ = value.(𝐩)

    xᵅ = solve(ZP, M, pᵥ; kwargs...)
    𝐱ᵅ = Dual{T}(xᵅ, one(xᵅ))

    fₓ = partials(f(𝐱ᵅ, pᵥ), 1)
    fₚ = partials(f(xᵅ, 𝐩))

    Dual{T}(xᵅ, - fₚ / fₓ)
end

end
