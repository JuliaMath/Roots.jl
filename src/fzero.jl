## MATLAB interfcae to find_zero
## Main functions are
## * fzero(f, ...) to find _a_ zero of f
## * fzeros(f, ...) to attempt to find all zeros of f



"""

Find zero of a function using an iterative algorithm

* `f`: a scalar function or callable object
* `x0`: an initial guess, finite valued.

Keyword arguments:

* `ftol`: tolerance for a guess `abs(f(x)) < ftol + max(1, |xnew|) * ftolrel`
* `ftolrel`: relative tolerance for convergence towards 0 of f(x)
* `xtol`: stop if `abs(xold - xnew) <= xtol + max(1, |xnew|)*xtolrel`
* `xtolrel`: see `xtol`
* `maxeval`: maximum number of steps
* `verbose`: Boolean. Set `true` to trace algorithm
* `order`: Can specify order of algorithm. 0 is most robust, also 1, 2, 5, 8, 16.
* `kwargs...` passed on to different algorithms. There are `maxfneval` when `order` is 1,2,5,8, or 16 and `beta` for orders 2,5,8,16,

This is a polyalgorithm redirecting to different algorithms based on the value of `order`. 

(The tolerance arguments can also be given through `atol`, `rtol`, `xatol`, and `xrtol`, as is done with `find_zero`.)

"""
function fzero(f, x0::Number; kwargs...)
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
fzero(f, a::Number, b::Number; kwargs...) = find_zero(f, (a,b), Bisection(); kwargs...)
fzero(f, bracket::Vector{T}; kwargs...)  where {T <: Number} = find_zero(f, bracket, Bisection(); kwargs...)
fzero(f, bracket::Tuple{T,S}; kwargs...)  where {T <: Number, S<:Number} = find_zero(f, bracket, Bisection();kwargs...)




"""
Find zero using Newton's method.
"""
fzero(f::Function, fp::Function, x0::Real; kwargs...) = newton(f, fp, float(x0); kwargs...)






# match fzero up with find_zero
@noinline function derivative_free(f, x0; order::Int=0, kwargs...) 
    
    if order == 0
        method = Order0()
    elseif order == 1
        method = Order1()
    elseif order == 2
        method = Order2()
    elseif order == 5
        method = Order5()
    elseif order == 8
        method = Order8()
    elseif order == 16
        method = Order16()
    else
        warn("Invalid order. Valid orders are 0, 1, 2, 5, 8, and 16")
        throw(ArgumentError())
    end

    d = Dict(kwargs)
    for (o, n) in ((:ftol, :atol), (:ftolrel, :rtol),
         (:xtol, :xatol), (:xtolrel, :xrtol))
        if haskey(d, o)
            d[n] = d[o]
        end
    end


    
    find_zero(f, x0, method; d...)
end

## 
"""
Find a zero within a bracket with an initial guess to *possibly* speed things along.
"""
function fzero(f, x0::Real, bracket::Vector{T}; kwargs...)  where {T <: Number}

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
function fzeros(f, a::Number, b::Number; kwargs...)  
    find_zeros(f, float(a), float(b); kwargs...)
end
fzeros(f, bracket::Vector{T}; kwargs...) where {T <: Number} = fzeros(f, a, b; kwargs...)
fzeros(f, bracket::Tuple{T,S}; kwargs...) where {T <: Number, S<:Number} = fzeros(f, a, b; kwargs...) 

