## MATLAB interfcae to find_zero
## Main functions are
## * fzero(f, ...) to find _a_ zero of f
## * fzeros(f, ...) to attempt to find all zeros of f



"""
    fzero(f, x0; order=0; kwargs...)

Find zero of a function using an iterative algorithm

* `f`: a scalar function or callable object
* `x0`: an initial guess, finite valued.
* `order`
    
This is a polyalgorithm redirecting to different algorithms based on the value of `order`. Dispatches to `find_zero(f, x0, OrderN(); kwargs...)`.

"""
function fzero(f, x0::Number; kwargs...)
    x = float(x0)
    isinf(x) && throw(ConvergenceFailed("An initial value must be finite"))
    derivative_free(f, x; kwargs...)
end



"""
    fzero(f, a, b; kwargs...)
    
Find zero of a function within a bracket, [a,b].

Dispatches to `find_zero(f, (a,b), Bisection())`.

"""
fzero(f, a::Number, b::Number; kwargs...) = find_zero(f, (a,b), Bisection(); kwargs...)
fzero(f, bracket::Vector{T}; kwargs...)  where {T <: Number} = find_zero(f, bracket, Bisection(); kwargs...)
fzero(f, bracket::Tuple{T,S}; kwargs...)  where {T <: Number, S<:Number} = find_zero(f, bracket, Bisection();kwargs...)




"""
    fzero(f, fp, x0; kwargs)
    
Find zero using Newton's method.

Dispatches to `find_zero((f,fp), x0, Roots.Newton(); kwargs...)`.
    
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


# ## 
# """
#     fzero(f, x0, bracket; kwargs...)
#    
# Find a zero within a bracket with an initial guess to *possibly* speed things along.
#
# Dispatches to the `A42` method.
#    
#"""
@deprecate fzero(f, x0::Real, bracket::Vector; kwargs...)  fzero(f, bracket)




## fzeros
"""

`fzeros(f, a, b; kwargs...)`
    
Searches for all zeros of `f` within an interval `(a,b)`. Assume neither `a` or `b` is a zero.

Dispatches to `find_zeros(f, a, b; kwargs...)`.
"""        
function fzeros(f, a::Number, b::Number; kwargs...)  
    find_zeros(f, float(a), float(b); kwargs...)
end
fzeros(f, bracket::Vector{T}; kwargs...) where {T <: Number} = fzeros(f, a, b; kwargs...)
fzeros(f, bracket::Tuple{T,S}; kwargs...) where {T <: Number, S<:Number} = fzeros(f, a, b; kwargs...) 

