module Roots

using Polynomials
import Polynomials: roots
export roots

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
include("real_roots.jl")

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
fzeros(p::Poly) = real_roots(p)

function fzero(p::Poly, x0::Real; kwargs...) 
    derivative_free(convert(Function, p), x0; kwargs...)
end
function fzero{T <: Real}(p::Poly, bracket::Vector{T}; kwargs...) 
    find_zero(convert(Function, p), bracket[1], bracket[2]; kwargs...)
end
function fzero{T <: Real}(p::Poly, x0::Real, bracket::Vector{T}; kwargs...) 
    fzero(convert(Function,p), x0, bracket; kwargs...)
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
    a,b = sort(map(float,bracket))
    try
        ex = secant_method_bracket(f, a, b, float(x0))
    catch ex
        if isa(ex, StateConverged) 
            return(ex.x0)
        else
            throw(ex)
        end
    end
end


## find all *simple* zeros in a bracket
function fzeros{T <: Real}(f::Function, bracket::Vector{T}; kwargs...) 
    ## check if a poly
    try
        filter(x -> bracket[1] <= x <= bracket[2], real_roots(convert(Poly, f)))
    catch e
        find_zeros(f, bracket[1], bracket[2]; kwargs...)
    end
end
fzeros(f::Function, a::Real, b::Real; kwargs...) = find_zeros(f, a, b; kwargs...)
function fzeros(f::Function)
    p = poly([0.0])
    try
        p = convert(Poly, f)
    catch e
        error("If f(x) is not a polynomial in x, then an interval to search over is needed")
    end
    real_roots(p)
end

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

