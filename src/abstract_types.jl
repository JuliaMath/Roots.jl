### Method types
abstract type AbstractUnivariateZeroMethod end
Base.broadcastable(method::AbstractUnivariateZeroMethod) = Ref(method)

abstract type AbstractBracketingMethod <: AbstractUnivariateZeroMethod end
abstract type AbstractBisectionMethod <: AbstractBracketingMethod end

abstract type AbstractNonBracketingMethod <: AbstractUnivariateZeroMethod end
abstract type AbstractSecantMethod <: AbstractNonBracketingMethod end

abstract type AbstractDerivativeMethod <: AbstractNonBracketingMethod end
abstract type AbstractNewtonLikeMethod <: AbstractDerivativeMethod end
abstract type AbstractHalleyLikeMethod <: AbstractDerivativeMethod end
abstract type AbstractΔMethod <: AbstractHalleyLikeMethod end

# deprecated but not clear way to do so, hence these definitions not to be used
const AbstractBracketing = AbstractBracketingMethod
const AbstractBisection = AbstractBisectionMethod
const AbstractNonBracketing = AbstractNonBracketingMethod
const AbstractSecant = AbstractSecantMethod

### State
abstract type AbstractUnivariateZeroState{T,S} end
