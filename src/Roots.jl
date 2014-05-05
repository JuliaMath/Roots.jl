module Roots
using Polynomial

export fzero, @fzero,
       fzeros, @fzeros,
       find_zero, 
       newton, halley,
       secant_method, steffensen,
       multroot,
       D, D2


## load in files
include("fzero.jl")
include("adiff.jl")
include("newton.jl")
include("derivative_free.jl")
include("multroot.jl")
include("SOLVE.jl")

## Main interface
##
## Dispatch handles many cases with reasonable algorithm
##
## args:
## f: Function of R -> R. May also be of `Poly` type.
##
## x0: initial guess for iterative algorithms. Required for non-poly, non-bracketed problems
##
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
function fzero(p::Poly, x0::Real; kwargs...) 
    derivative_free(convert(Function, p), x0; kwargs...)
end
function fzero{T <: Real}(p::Poly, bracket::Vector{T}; kwargs...) 
    find_zero(convert(Function, p), bracket[1], bracket[2]; kwargs...)
end
function fzero{T <: Real}(p::Poly, x0::Real, bracket::Vector{T}; kwargs...) 
    derivative_free_bracket(convert(Function,p), x0, bracket; kwargs...)
end

##
newton(p::Poly, x0::Real; kwargs...) = newton(convert(Function, p), convert(Function, polyder(p)), x0; kwargs...)
##

## Functions

## bracket
function fzero(f::Function, a::Real, b::Real; kwargs...) 
    find_zero(f, promote(float(a), b)...; kwargs...)
end
function fzero{T <: Real}(f::Function, bracket::Vector{T}; kwargs...) 
    fzero(f, bracket[1], bracket[2]; kwargs...)
end
function fzero{T <: Real}(f::Function, x0::Real, bracket::Vector{T}; kwargs...) 
    derivative_free_bracket(f, x0, bracket; kwargs...)
end

## simplify calling of function (no x ->, rather @fzero...)
macro fzero(expr::Expr, a, b)
    quote
        fzero(x -> $(expr), $(a), $(b))
    end
end

macro fzero(expr::Expr, a)
    quote
        fzero(x -> $(expr), $(a))
    end
end


## find all zeros in a bracket
function fzeros{T <: Real}(f::Function, bracket::Vector{T}; kwargs...) 
    find_zeros(f, bracket[1], bracket[2]; kwargs...)
end
fzeros(f::Function, a::Real, b::Real; kwargs...) = find_zeros(f, a, b; kwargs...)


macro fzeros(expr::Expr, a, b)
    quote
        fzeros(x -> $(expr), $(a), $(b))
    end
end


## derivative free
fzero(f::Function, x0::Real; kwargs...) = derivative_free(f, x0; kwargs...)

## use newton if two functions passed
fzero(f::Function, fp::Function, x0::Real; kwargs...) = newton(f, fp, x0; kwargs...)

## use Automatic differentiation here
newton(f::Function, x::Real; kwargs...) =  newton(f, D(f), x; kwargs...)
halley(f::Function, x::Real; kwargs...) = halley(f, D(f), D2(f), x; kwargs...)


end

