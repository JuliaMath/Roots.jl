### Method types
abstract type AbstractUnivariateZeroMethod end

abstract type AbstractBracketingMethod <: AbstractUnivariateZeroMethod end
abstract type AbstractBisectionMethod <: AbstractBracketingMethod end

abstract type AbstractNonBracketingMethod <: AbstractUnivariateZeroMethod end
abstract type AbstractSecantMethod <: AbstractNonBracketingMethod end

abstract type AbstractDerivativeMethod <:AbstractNonBracketingMethod  end
abstract type AbstractNewtonLikeMethod <: AbstractDerivativeMethod end
abstract type AbstractHalleyLikeMethod <: AbstractDerivativeMethod  end
abstract type AbstractΔMethod <: AbstractHalleyLikeMethod end



### State
abstract type AbstractUnivariateZeroState{T,S} end
