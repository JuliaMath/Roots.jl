using Roots
using Test
import SpecialFunctions.erf

include("./test_find_zero.jl")
include("./test_bracketing.jl")
include("./test_derivative_free.jl")
include("./test_simple.jl")
include("./test_find_zeros.jl")
include("./test_fzero.jl")
include("./test_newton.jl")
include("./test_simple.jl")

#include("./test_composable.jl")
#include("./runbenchmarks.jl")
#include("./test_derivative_free_interactive.jl")
