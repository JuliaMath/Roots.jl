using Base.Test
import Roots.newton, Roots.halley

@test_approx_eq_eps newton(sin, cos, 0.5) 0.0 100*eps(1.0)
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


## tests with auto derivatives
@test_approx_eq newton(sin, 3) float(pi)
@test_approx_eq halley(sin, 3) float(pi)

## More tests with autoderivaitves. Derivative of D operation:
isdefined(ForwardDiff, :derivative) && @test_approx_eq newton(D(sin), 1.5) pi/2

## test with Complex input

@test_approx_eq real(newton(x ->  x^3 - 1, x ->  3x^2, 1+im)) 1
@test_approx_eq real(newton(x ->  x^3 - 1, x ->  3x^2, 1+10im)) (-1/2)
