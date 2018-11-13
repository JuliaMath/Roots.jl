## MATLAB interfcae to find_zero
## Main functions are
## * fzero(f, ...) to find _a_ zero of f, a univariate function
## * fzeros(f, ...) to attempt to find all zeros of f, a univariate function
## unlike `find_zero` these do not specialize on f, so
## will be faster the first use, and slower for subsequent uses (for the same f)

"""
    fzero(f, x0; order=0; kwargs...)

Find zero of a function using an iterative algorithm

* `f`: a scalar function or callable object
* `x0`: an initial guess, finite valued.
* `order`: An integer, symbol, or string indicating the algorithm to
   use. The default is `Order0`.

"""
function fzero(f, x0::Number; kwargs...)
    x = float(x0)
    isinf(x) && throw(ConvergenceFailed("An initial value must be finite"))
    derivative_free(f, x; kwargs...)
end


"""
    fzero(f, x0, M, [N]; kwargs...)
"""
function fzero(f, x0, M::AbstractUnivariateZeroMethod; kwargs...)
    find_zero(FnWrapper(f), x0, M; kwargs...)
end

function fzero(f, x0, M::AbstractUnivariateZeroMethod, N::AbstractBracketing; kwargs...)
    find_zero(FnWrapper(f), x0, M, N; kwargs...)
end

"""
    fzero(f, a, b; kwargs...)

If `order` is not specified, will use bisection to find a zero of a
function within a bracket, [a,b]. Otherwise, will use `x0=(a,b)` for
the method specified by `order`.

"""
function fzero(f, bracket::Tuple{T,S}; kwargs...)  where {T <: Number, S<:Number}
    if haskey(kwargs, :order)
        find_zero(FnWrapper(f), bracket, _method_lookup(kwargs,:order); kwargs...)
    else
        find_zero(FnWrapper(f), bracket, Bisection();kwargs...)
    end
end

fzero(f, a::Number, b::Number; kwargs...) = fzero(f, (a,b); kwargs...)
fzero(f, bracket::Vector{T}; kwargs...)  where {T <: Number} = fzero(f,(bracket[1],bracket[2]); kwargs...)




"""
    fzero(f, fp, x0; kwargs)

Find zero using Newton's method.

Dispatches to `find_zero((f,fp), x0, Roots.Newton(); kwargs...)`.

"""
fzero(f::Function, fp::Function, x0::Real; kwargs...) = find_zero((f,fp), x0, Newton(); kwargs...)





# match fzero up with find_zero
_method_lookup = Dict(0   => Order0(),
                     :0  => Order0(),
                     "0" => Order0(),

                      1    => Order1(),
                      :1   => Order1(),
                      "1"  => Order1(),
                      :secant => Order1(),
                      "1B" => Order1B(),
                      :king => Order1B(),

                     2    => Order2(),
                     :2   => Order2(),
                     :steffensen   => Order2(),
                     "2"  => Order2(),
                     "2B" => Order2B(),
                     :esser   => Order2B(),

                      5  => Order5(),
                     :5  => Order5(),
                     "5" => Order5(),

                     8   => Order8(),
                     :8  => Order8(),
                     "8" => Order8(),

                     16   => Order16(),
                     :16  => Order16(),
                     "16" => Order16(),
)

@noinline function derivative_free(f, x0; order=0, kwargs...)


    if haskey(_method_lookup, order)
        M = _method_lookup[order]
    else
        warn("Invalid order. Valid orders are 0, 1, 2, 5, 8, and 16")
        throw(ArgumentError())
    end

    # d = (kv[1] == :ftol ? :atol=>kv[2] :
    #      kv[1] == :ftolrel ? :rtol=>kv[2] :
    #      kv[1] == :xtol ? :xatol=>kv[2] :
    #      kv[1] == :xtolrel ? xrtol=>kv[2] :
    #      kv[1] => kv[1] for kv in kwargs)

    d = Dict(kwargs)
     for (o, n) in ((:ftol, :atol), (:ftolrel, :rtol),
                    (:xtol, :xatol), (:xtolrel, :xrtol))
         if haskey(d, o)
             d[n] = d[o]
         end
     end

    find_zero(FnWrapper(f), x0, M; d...)
end






## fzeros
"""

`fzeros(f, a, b; kwargs...)`

Searches for all zeros of `f` within an interval `(a,b)`. Assume neither `a` or `b` is a zero.

Dispatches to `find_zeros(f, a, b; kwargs...)`.
"""
function fzeros(f, a::Number, b::Number; kwargs...)
    find_zeros(FnWrapper(f), float(a), float(b); kwargs...)
end
fzeros(f, bracket::Vector{T}; kwargs...) where {T <: Number} = fzeros(f, a, b; kwargs...)
fzeros(f, bracket::Tuple{T,S}; kwargs...) where {T <: Number, S<:Number} = fzeros(f, a, b; kwargs...)
