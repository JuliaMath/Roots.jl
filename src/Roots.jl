module Roots

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 1
end

using Printf

export fzero,
       fzeros,
       secant_method

export find_zero, find_zero!, find_zeros,
       Order0, Order1, Order2, Order5, Order8, Order16

export Bisection, FalsePosition

## load in files
include("utils.jl")
include("find_zero.jl")
include("bracketing.jl")
include("derivative_free.jl")
include("simple.jl")
include("find_zeros.jl")
include("newton.jl")
include("fzero.jl")



end
