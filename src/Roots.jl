module Roots

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 1
end

using Printf
import CommonSolve
import CommonSolve: solve, solve!, init
using Setfield
import ChainRulesCore

export fzero, fzeros, secant_method

export find_zero,
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

include("utils.jl")
include("abstract_types.jl")
include("state.jl")
include("convergence.jl")
include("functions.jl")
include("trace.jl")
include("find_zero.jl")
include("hybrid.jl")
include("chain_rules.jl")

include("Bracketing/bracketing.jl")
include("Bracketing/bisection.jl")
include("Bracketing/alefeld_potra_shi.jl")
include("Bracketing/brent.jl")
include("Bracketing/ridders.jl")
include("Bracketing/itp.jl")
include("Bracketing/chandrapatlu.jl")
include("Bracketing/false_position.jl")

include("DerivativeFree/derivative_free.jl")
include("DerivativeFree/secant.jl")
include("DerivativeFree/steffensen.jl")
include("DerivativeFree/order5.jl")
include("DerivativeFree/order8.jl")
include("DerivativeFree/order16.jl")
include("DerivativeFree/king.jl")
include("DerivativeFree/esser.jl")
include("DerivativeFree/order0.jl")

include("Derivative/newton.jl")
include("Derivative/halley_like.jl")
include("Derivative/thukralb.jl")
include("Derivative/lith.jl")

include("find_zeros.jl")
include("simple.jl")
include("alternative_interfaces.jl")

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
