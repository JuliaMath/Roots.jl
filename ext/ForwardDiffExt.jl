module ForwardDiffExt

using Roots
using ForwardDiff
using ForwardDiffChainRules

@ForwardDiff_frule solve(ZP, M, p::T) where {T <: ForwardDiff.Dual}
@ForwardDiff_frule solve(ZP, M, p::AbstractVector{T}) where {T <: ForwardDiff.Dual}





end
