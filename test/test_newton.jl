using Test
import Roots.newton, Roots.halley

@test_approx_eq newton(sin, cos, 0.5) 0.0
@test_approx_eq newton(cos, x -> -sin(x), 1.0) pi/2
@test_approx_eq newton(x -> x^2 - 2x - 1, x -> 2x - 2, 3.0) 2.414213562373095
@test_approx_eq_eps newton(x -> exp(x) - cos(x),
                           x -> exp(x) + sin(x), 3.0) 0.0 eps(1.0) 
@test_approx_eq halley(x -> x^2 - 2x - 1,
                       x -> 2x - 2,
                       x -> 2, 3.0) 2.414213562373095
@test_approx_eq_eps halley(x -> exp(x) - cos(x),
                           x -> exp(x) + sin(x),
                           x -> exp(x) + cos(x), 3.0) 0.0 eps(1.0)
