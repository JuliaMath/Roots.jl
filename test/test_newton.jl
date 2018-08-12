using Test

@testset "Newton, halley" begin
    
    @test abs(Roots.newton(sin, cos, 0.5) - 0.0) <= 100*eps(1.0)
    @test Roots.newton(cos, x -> -sin(x), 1.0)  ≈ pi/2 
    @test Roots.newton(x -> x^2 - 2x - 1, x -> 2x - 2, 3.0)  ≈ 2.414213562373095
    @test abs(Roots.newton(x -> exp(x) - cos(x), x -> exp(x) + sin(x), 3.0) - 0.0) <= 1e-14
    @test Roots.halley(x -> x^2 - 2x - 1,x -> 2x - 2,x -> 2, 3.0)  ≈ 2.414213562373095
    a = Roots.halley(x -> exp(x) - cos(x),
                     x -> exp(x) + sin(x),
                     x -> exp(x) + cos(x), 3.0)
    @test abs(a - 0.0) <= 1e-14
    
    
    
    ## test with Complex input
    
    @test real(Roots.newton(x ->  x^3 - 1, x ->  3x^2, 1+im)) ≈ 1.0
    @test real(Roots.newton(x ->  x^3 - 1, x ->  3x^2, 1+10im)) ≈ (-1/2)
end
