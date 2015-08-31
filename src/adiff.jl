## add in simple derivative operator using PowerSeries package
## D(f), D(f,k) returns the derivative to power k (k=1 is default)

"""

Take derivative of order `k` of a function.

Arguments:
* `f::Function`: a mathematical function from R to R.
* `k::Int=1`: A non-negative integer specifying the order of the
  derivative. Values larger than 8 can be slow to compute.

Uses the `PowerSeries` package for forward automatic
differentiation. This restricts which functions are possible to
differentiate. In particular, functions returned by `D` are not able
to be used with this operator. To find higher-order derivatives, pass
in values of `k` larger than 1.

"""
function D(f::Function, k::Int=1)
    k == 0 && return(f)
    k < 0 && error("The order of the derivative must be non-negative")

    if VERSION < v"0.4-"
        tmp = ntuple(k-1, z->0.0)
    else
        tmp = ntuple(z -> 0.0, k-1)
    end

    function(x)
        y = series(tuple(x, 1.0, tmp...)...)
        factorial(k) * f(y).(symbol(pop!(fieldnames(y))))
    end
end

D2(f::Function) = D(f, 2)
