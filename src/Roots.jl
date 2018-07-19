__precompile__(true)
module Roots


if VERSION >= v"0.7-"
    using Printf
else
    using Missings
end

using Compat: @nospecialize, lastindex, range


export fzero,
       fzeros
       #newton, halley,  # deprecate these 4?
       #secant_method, steffensen

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
