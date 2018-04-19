__precompile__(true)
module Roots


if VERSION >= v"0.7-"
    using Printf
else
    using Missings
end

using ForwardDiff
using Compat: @nospecialize


export fzero,
       fzeros,
       newton, halley,  # deprecate these 4?
       secant_method, steffensen,
       D

export find_zero,
       Order0, Order1, Order2, Order5, Order8, Order16
export Bisection, FalsePosition

## load in files
include("adiff.jl")
include("utils.jl")
include("find_zero.jl")
include("bracketing.jl")
include("derivative_free.jl")
include("newton.jl")


## Main functions are
## * fzero(f, ...) to find _a_ zero of f
## * fzeros(f, ...) to attempt to find all zeros of f



"""

Find zero of a function using an iterative algorithm

* `f`: a scalar function or callable object
* `x0`: an initial guess, finite valued.

Keyword arguments:

* `ftol`: tolerance for a guess `abs(f(x)) < ftol`
* `xtol`: stop if `abs(xold - xnew) <= xtol + max(1, |xnew|)*xtolrel`
* `xtolrel`: see `xtol`
* `maxeval`: maximum number of steps
* `verbose`: Boolean. Set `true` to trace algorithm
* `order`: Can specify order of algorithm. 0 is most robust, also 1, 2, 5, 8, 16.
* `kwargs...` passed on to different algorithms. There are `maxfneval` when `order` is 1,2,5,8, or 16 and `beta` for orders 2,5,8,16,

This is a polyalgorithm redirecting different algorithms based on the value of `order`. 

"""
function fzero(f, x0::Real; kwargs...)
    x = float(x0)
    isinf(x) && throw(ConvergenceFailed("An initial value must be finite"))
    derivative_free(f, x; kwargs...)
end



"""
Find zero of a function within a bracket

Uses a modified bisection method for non `big` arguments

Arguments:

* `f` A scalar function or callable object
* `a` left endpont of interval
* `b` right endpont of interval
* `xtol` optional additional tolerance on sub-bracket size.

For a bracket to be valid, it must be that `f(a)*f(b) < 0`.

For `Float64` values, the answer is such that `f(prevfloat(x)) *
f(nextfloat(x)) < 0` unless a non-zero value of `xtol` is specified in
which case, it stops when the sub-bracketing produces an bracket with
length less than `max(xtol, abs(x1)*xtolrel)`. 

For `Big` values, which defaults to the algorithm of Alefeld, Potra, and Shi, a
default tolerance is used for the sub-bracket length that can be enlarged with `xtol`.

If `a==-Inf` it is replaced with `nextfloat(a)`; if `b==Inf` it is replaced with `prevfloat(b)`.

Example:

    `fzero(sin, 3, 4)` # find pi
    `fzero(sin, [big(3), 4]) find pi with more digits
"""
function fzero(f, a::Real, b::Real; kwargs...)
    bracket = adjust_bracket((a, b))
    a0, b0 = a < b ? promote(float(a), b) : promote(float(b), a)
    (a0 == -Inf) && (a0 = nextfloat(a0))
    (b0 == Inf) && (b0 = prevfloat(b0))
    find_zero(f, bracket, Bisection(); kwargs...)
end

"""
Find a zero with bracket specified via `[a,b]`, as `fzero(sin, [3,4])`.
"""
function fzero(f, bracket::Vector{T}; kwargs...)  where {T <: Real}
    fzero(f, float(bracket[1]), float(bracket[2]); kwargs...)
end


"""
Find a zero within a bracket with an initial guess to *possibly* speed things along.
"""
function fzero(f, x0::Real, bracket::Vector{T}; kwargs...)  where {T <: Real}

    a, b = adjust_bracket(bracket)

    try
        ex = a42a(f, a, b, float(x0); kwargs...)
    catch ex
        if isa(ex, StateConverged) 
            return(ex.x0)
        else
            rethrow(ex)
        end
    end
end


"""
Find zero using Newton's method.
"""
fzero(f::Function, fp::Function, x0::Real; kwargs...) = newton(f, fp, float(x0); kwargs...)



## fzeros
"""

`fzeros(f, a, b)`
    
Attempt to find all zeros of `f` within an interval `[a,b]`.

Simple algorithm that splits `[a,b]` into subintervals and checks each
for a root.  For bracketing subintervals, bisection is
used. Otherwise, a derivative-free method is used. If there are a
large number of zeros found relative to the number of subintervals, the
number of subintervals is increased and the process is re-run.

There are possible issues with close-by zeros and zeros which do not
cross the origin (non-simple zeros). Answers should be confirmed
graphically, if possible.

"""        
function fzeros(f, a::Real, b::Real; kwargs...)  
    find_zeros(f, float(a), float(b); kwargs...)
end
fzeros(f, bracket::Vector{T}; kwargs...) where {T <: Real} = fzeros(f, a, b; kwargs...) 



end
