using Roots
using Compat.Test
import SpecialFunctions.erf

include("./test_find_zero.jl")
include("./test_fzero.jl")
include("./test_find_zeros.jl")
include("./test_newton.jl")
include("./test_simple.jl")
include("./RootTesting.jl")

#include("./test_composable.jl")

#run_benchmark_tests()

#include("./test_fzero3.jl")
#run_robustness_test()

#include("./test_derivative_free.jl")
