module Roots

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 1
end

using Printf
import CommonSolve
import CommonSolve: solve, solve!, init
using Setfield

export fzero, fzeros, secant_method

export find_zero,
    find_zero!,
    find_zeros,
    ZeroProblem,
    solve,
    solve!,
    init,
    Order0,
    Secant,
    Order1,
    Orderφ,
    Steffensen,
    Order2,
    Order5,
    Order8,
    Order16

export Bisection, A42, AlefeldPotraShi, FalsePosition

## load in files
include("utils.jl")
include("find_zero.jl")
include("bracketing.jl")
include("derivative_free.jl")
include("order0.jl")
include("simple.jl")
include("find_zeros.jl")
include("newton.jl")
include("lith.jl")
include("non_simple.jl")
include("fzero.jl")

# cf. https://github.com/JuliaDocs/Documenter.jl/pull/1664/files
function _update_module_doc()
    path = joinpath(@__DIR__, "..", "README.md")
    text = read(path, String)
    # The code blocks in the README.md should be julia blocks for the syntax highlighter.
    text = replace(text, "```julia" => "```jldoctest readme")
    @doc text Roots
end
_update_module_doc()

end
