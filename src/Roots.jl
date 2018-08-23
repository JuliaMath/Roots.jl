VERSION < v"0.7.0-beta2.199" && __precompile__()
module Roots


if VERSION >= v"0.7-"
    using Printf
end

using Compat: @nospecialize, lastindex, range, Nothing


export fzero,
       fzeros,
       secant_method

export find_zero, find_zeros,
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
