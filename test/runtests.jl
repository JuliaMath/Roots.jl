using Roots
using Test
using Aqua

import SpecialFunctions.erf

struct SomeInterval{T}
    a::T
    b::T
end
SomeInterval(a, b) = SomeInterval(promote(a, b)...)
Base.extrema(I::SomeInterval) = I.a < I.b ? (I.a, I.b) : (I.b, I.a)

# count function calls
mutable struct Cnt
    cnt::Int
    f
    Cnt(f) = new(0, f)
end
(f::Cnt)(x) = (f.cnt += 1; f.f(x))
avg(x) = sum(x) / length(x)

include("./test_find_zero.jl")
include("./test_bracketing.jl")
include("./test_derivative_free.jl")
include("./test_simple.jl")
include("./test_find_zeros.jl")
include("./test_fzero.jl")
include("./test_newton.jl")
include("./test_chain_rules.jl")
include("./test_simple.jl")

include("./test_composable.jl")
VERSION >= v"1.6.0" && include("./test_allocations.jl")

#include("./runbenchmarks.jl")
#include("./test_derivative_free_interactive.jl")

Aqua.test_all(Roots; ambiguities=false, project_toml_formatting=false)
