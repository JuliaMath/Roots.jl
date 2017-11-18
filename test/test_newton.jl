using Compat.Test
import Roots.newton, Roots.halley

@test norm(newton(sin, cos, 0.5) - 0.0) <= 100*eps(1.0)
@test newton(cos, x -> -sin(x), 1.0)  ≈ pi/2 
@test newton(x -> x^2 - 2x - 1, x -> 2x - 2, 3.0)  ≈ 2.414213562373095
@test norm(newton(x -> exp(x) - cos(x), x -> exp(x) + sin(x), 3.0) - 0.0) <= 1e-14
@test halley(x -> x^2 - 2x - 1,x -> 2x - 2,x -> 2, 3.0)  ≈ 2.414213562373095
a = halley(x -> exp(x) - cos(x),
           x -> exp(x) + sin(x),
           x -> exp(x) + cos(x), 3.0)
@test norm(a - 0.0) <= 1e-14


## tests with auto derivatives
@test newton(sin, 3) ≈ pi
@test halley(sin, 3) ≈ pi

## More tests with autoderivaitves. Derivative of D operation:
isdefined(ForwardDiff, :derivative) && @test newton(D(sin), 1.5) ≈ pi/2

## test with Complex input

@test real(newton(x ->  x^3 - 1, x ->  3x^2, 1+im)) ≈ 1.0
@test real(newton(x ->  x^3 - 1, x ->  3x^2, 1+10im)) ≈ (-1/2)
