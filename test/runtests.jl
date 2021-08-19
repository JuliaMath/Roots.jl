using Roots
using Test
import SpecialFunctions.erf

struct SomeInterval
    a
    b
end
Base.extrema(I::SomeInterval) = I.a < I.b ? (I.a, I.b) : (I.b, I.a)


include("./test_find_zero.jl")
include("./test_bracketing.jl")
include("./test_derivative_free.jl")
include("./test_simple.jl")
include("./test_find_zeros.jl")
include("./test_fzero.jl")
include("./test_newton.jl")
include("./test_simple.jl")

include("./test_composable.jl")
include("./test_allocations.jl")

#include("./runbenchmarks.jl")
#include("./test_derivative_free_interactive.jl")
