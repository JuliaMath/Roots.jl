__precompile__(true)
module Roots


import Base: *


#using ForwardDiff
#using Compat

export fzero,
       fzeros,
#       newton, halley,
       secant_method, steffensen
#       D
#export multroot, D2  # deprecated

export find_zero,
       Order0, Order1, Order2, Order5, Order8, Order16
export Bisection, FalsePosition
# export Bisection, Secant, Steffensen, Newton, Halley

## load in files
#include("adiff.jl")
include("find_zero.jl")
include("bracketing.jl")
include("derivative_free.jl")
#include("newton.jl")


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
    a,b = sort([a,b])
    a,b = promote(float(a), b)
    (a == -Inf) && (a = nextfloat(a))
    (b == Inf) && (b = prevfloat(b))
    (isinf(a) | isinf(b)) && throw(ConvergenceFailed("A bracketing interval must be bounded"))
    find_zero(f, [a,b], Bisection(); kwargs...)
end

"""
Find a zero with bracket specified via `[a,b]`, as `fzero(sin, [3,4])`.
"""
function fzero{T <: Real}(f, bracket::Vector{T}; kwargs...) 
    fzero(f, float(bracket[1]), float(bracket[2]); kwargs...)
end

"""
Find a zero within a bracket with an initial guess to *possibly* speed things along.
"""
function fzero{T <: Real}(f, x0::Real, bracket::Vector{T}; kwargs...) 
    a,b = sort(map(float,bracket))
    (a == -Inf) && (a = nextfloat(a))
    (b == Inf) && (b = prevfloat(b))
    (isinf(a) | isinf(b)) && throw(ConvergenceFailed("A bracketing interval must be bounded"))    
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
fzeros{T <: Real}(f, bracket::Vector{T}; kwargs...)  = fzeros(f, a, b; kwargs...)


## deprecate Polynomial calls

## Don't want to load `Polynomials` to do this...
# @deprecate roots(p) fzeros(p, a, b) 

@deprecate D2(f) D(f,2)
fzeros(p) = Base.depwarn("""
Calling fzeros with just a polynomial is deprecated.
Either:
   * Specify an interval to search over: fzeros(p, a, b).
   * Use the `poly_roots(p, Over.R)` call from `PolynomialZeros`                         
   * Use `Polynomials` or `PolynomialRoots` and filter. For example, 

```
using Polynomials
x=variable()
filter(isreal, roots(f(x)))
```
                         
""",:fzeros)

multroot(p) = Base.depwarn("The multroot function has moved to the PolynomialZeros package.",:multroot)

polyfactor(p) = Base.depwarn("""
The polyfactor function is in the PolynomialFactors package:

* if `p` is of type `Poly`: `PolynomialFactors.factor(p)`.
                             
* if `p` is a function, try:
```
using PolynomialFactors, Polynomials
x = variable(Int)
PolynomialFactors.factor(p(x))
```
                             
""",:polyfactor)

end
