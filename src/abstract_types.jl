### Method types
abstract type AbstractUnivariateZeroMethod end

abstract type AbstractBracketing <: AbstractUnivariateZeroMethod end
abstract type AbstractBisection <: AbstractBracketing end
abstract type AbstractAcceleratedBisection <: AbstractBisection end

abstract type AbstractNonBracketing <: AbstractUnivariateZeroMethod end
abstract type AbstractSecant <: AbstractNonBracketing end

abstract type AbstractNewtonLikeMethod <: AbstractNonBracketing  end
abstract type AbstractHalleyLikeMethod <: AbstractNonBracketing  end
abstract type AbstractΔMethod <: AbstractHalleyLikeMethod end



### State
abstract type AbstractUnivariateZeroState{T,S} end
