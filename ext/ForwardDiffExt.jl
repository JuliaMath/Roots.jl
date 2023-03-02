module ForwardDiffExt

using Roots
using ForwardDiff
using ForwardDiffChainRules

@ForwardDiff_frule solve(ZP::ZeroProblem, M::Roots.AbstractUnivariateZeroMethod,
                         p::ForwardDiff.Dual)
@ForwardDiff_frule solve(ZP::ZeroProblem, M::Roots.AbstractUnivariateZeroMethod,
                         p::AbstractVector{<: ForwardDiff.Dual})

end
