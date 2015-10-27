## add in simple derivative operator using the ForwardDiff package
## D(f), D(f,k) returns the derivative to power k (k=1 is default)

"""

Take derivative of order `k` of a function.

Arguments:
* `f::Function`: a mathematical function from R to R.
* `k::Int=1`: A non-negative integer specifying the order of the
  derivative. Values larger than 8 can be slow to compute.

Wrapper around `derivative` function in `ForwardDiff`

"""
function D(f::Function, k::Int=1)
    k < 0 && error("The order of the derivative must be non-negative")
    k == 0 && return(f)
    k == 1 && return(derivative(f))
    D(derivative(f))
end

D2(f::Function) = D(f, 2)
