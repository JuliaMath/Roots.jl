module Roots
using Polynomial

export fzero,
       find_zero, newton, halley,
       thukral,
       multroot,
       D, D2


## load in files
include("fzero.jl")
include("adiff.jl")
include("newton.jl")
include("thukral.jl")
include("multroot.jl")

## Main interface
##
## Dispatch handles many cases with most rapidly converging algorithm
##
## args:
## f: Function of R -> R. May also be of `Poly` type.
## x0: initial guess for iterative algorithms. Required for non-poly, non-bracketed problems
## bracket: bracket [a,b] with f(a) * f(b) < 0. Bracketing guarantees a root will be found in the interval.
##
## kwargs: 
##
## tol: test if | f(x_n) | <= tol
## delta: test if | x_{n+1} - x_n } <= delta
## max_iter: maximum number of steps before giving up on convergene
## 
## Poly types
fzero(p::Poly) = multroot(p)
fzero(p::Poly, x0::Real; kwargs...) = thukral(convert(Function, p), x0; kwargs...)
fzero{T <: Real}(p::Poly, bracket::Vector{T}; kwargs...) = find_zero(convert(Function, p), bracket[1], bracket[2]; kwargs...)
fzero{T <: Real}(p::Poly, x0::Real, bracket::Vector{T}; kwargs...) = thukral_bracket(convert(Function,p), x0, bracket; kwargs...)
##
newton(p::Poly, x0::Real; kwargs...) = newton(convert(Function, p), convert(Function, polyder(p)), x0; kwargs...)
##
## Functions
fzero(f::Function, x0::Real; kwargs...) = thukral(f, x0; kwargs...)
fzero(f::Function, a::Real, b::Real; kwargs...) = find_zero(f, a, b; kwargs...)
fzero{T <: Real}(f::Function, bracket::Vector{T}; kwargs...) = find_zero(f, bracket[1], bracket[2]; kwargs...)
fzero{T <: Real}(f::Function, x0::Real, bracket::Vector{T}; kwargs...) = thukral_bracket(f, x0, bracket; kwargs...)

## use newton if two functions passed
fzero(f::Function, fp::Function, x0::Real; kwargs...) = newton(f, fp, x0; kwargs...)


newton(f::Function, x::Real; kwargs...) =  newton(f, D(f), x; kwargs...)
halley(f::Function, x::Real; kwargs...) = halley(f, D(f), D2(f), x; kwargs...)


end

