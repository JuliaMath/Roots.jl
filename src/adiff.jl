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
    k == 0 && return(x -> f(float(x)))

     if isdefined(ForwardDiff, :derivative)
            D(derivative(f), k-1)
     else
         if k == 1
             g = forwarddiff_gradient(v->f(v[1]), Float64)
             x -> g([float(x)])[1]
         elseif k==2
             g1(x) = f(x[1])
             g2 = forwarddiff_hessian(g1, Float64)
             x -> g2([float(x), 0.0])[1,1]
         else
             error("Automatic derivatives for this version of ForwardDiff.jl can only handle k=0,1, or 2")
         end
     end
end

D2(f::Function) = D(f, 2)

## This conflicts with a definition in Calculus, but is better. 
## Base.ctranspose(f::Function) = D(f)
