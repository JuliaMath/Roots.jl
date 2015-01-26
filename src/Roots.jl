module Roots

using Polynomials
import Polynomials: roots

using PowerSeries

## * Docile is used for documentation
if VERSION < v"0.4.0-dev"
    using Docile
end
@document

export roots

export fzero,
       fzeros,
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


## Main functions are
## * fzero(f, ...) to find _a_ zero of f
## * fzeros(f, ...) to attempt to find all zeros of f



"""

Find zero of a function

* `f`: a scalar function 
* `x0`: an initial guess

Keyword arguments:

* `tol`: tolerance for a guess `abs(f(x)) < tol`
* `reltol`: relative tolerance
* `delta`: stop if `abs(xold - xnew) <= delta`
* `max_iter`: maximum number of steps
* `verbose`: Boolean. Set `true` to trace algorithm
* `order`: Can specify order of algorithm
* `kwargs...` passed on to different algorithms

"""
fzero(f::Function, x0::Real; kwargs...) = derivative_free(f, float(x0); kwargs...)


"""
Find zero of a function within a bracket

Uses a modified bisection method for non `big` arguments

Arguments:

* `f` function
* `a` left endpont of interval
* `b` right endpont of interval

For a bracket to be valid, it must be that `f(a)*f(b) < 0`.
"""
function fzero(f::Function, a::Real, b::Real; kwargs...) 
    find_zero(f, promote(float(a), b)...; kwargs...)
end

"""
Find a zero with bracket specified via `[a,b]`
"""
function fzero{T <: Real}(f::Function, bracket::Vector{T}; kwargs...) 
    fzero(f, bracket[1], bracket[2]; kwargs...)
end

"""
Find a zero within a bracket with an initial guess to *possibly* speed things along.
"""
function fzero{T <: Real}(f::Function, x0::Real, bracket::Vector{T}; kwargs...) 
    a,b = sort(map(float,bracket))
    try
        ex = secant_method_bracket(f, a, b, float(x0))
    catch ex
        if isa(ex, StateConverged) 
            return(ex.x0)
        else
            rethrow(ex)
        end
    end
end


"""

Find zero using Newton's method

"""
fzero(f::Function, fp::Function, x0::Real; kwargs...) = newton(f, fp, x0; kwargs...)



"""

Find zero of a polynomial with derivative free algorithm.

Arguments:

* `p` a `Poly` object
* `x0` an initial guess

Returns:

An approximate root or an error.

"""
function fzero(p::Poly, x0::Real; kwargs...) 
    derivative_free(convert(Function, p), x0; kwargs...)
end

function fzero{T <: Real}(p::Poly, bracket::Vector{T}; kwargs...) 
    find_zero(convert(Function, p), bracket[1], bracket[2]; kwargs...)
end
function fzero{T <: Real}(p::Poly, x0::Real, bracket::Vector{T}; kwargs...) 
    fzero(convert(Function,p), x0, bracket; kwargs...)
end


## fzeros

"""

Find real zeros of a polynomial

args:

`f`: Function of R -> R. May also be of `Poly` type.

`x0`: initial guess for iterative algorithms. Required for non-poly, non-bracketed problems

bracket: bracket [a,b] with f(a) * f(b) < 0. Bracketing guarantees a root will be found in the interval.

kwargs: 

`tol`: test if | f(x_n) | <= tol
`delta`: test if | x_{n+1} - x_n } <= delta
`max_iter`: maximum number of steps before giving up on convergene


"""
fzeros(p::Poly) = real_roots(p)


"""
Attempt to find all simple zeroes of `f` within an interval `[a,b]`.

Simple algorithm simply splits `[a,b]` into 250 subintervals and checks each for a bracket.
"""        
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






## use Automatic differentiation here
newton(f::Function, x::Real; kwargs...) =  newton(f, D(f), x; kwargs...)

##
newton(p::Poly, x0::Real; kwargs...) = newton(convert(Function, p), convert(Function, polyder(p)), x0; kwargs...)
##

halley(f::Function, x::Real; kwargs...) = halley(f, D(f), D2(f), x; kwargs...)


end

