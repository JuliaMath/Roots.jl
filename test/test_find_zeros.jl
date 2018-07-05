using Roots
using Compat.Test


@testset "find_zeros" begin

    rts = find_zeros(x -> exp(x) - x^4, -5, 20)
    @test length(rts) == 3

    n = 10
    rts = find_zeros(x -> cos(x)/n, 0, 10)
    @test length(rts) == 3

    W(n) = x -> prod(x-i for i in 1:n)
    rts = find_zeros(W(20), -1, 21)
    @test length(rts) == 20

    W1(n) = x -> prod((x-i)^i for i in 1:n)
    rts = find_zeros(W1(5), -0.5, 5.5, rtol=1e-12) # o/w even powers hard
    @test length(rts) == 5

    T11(x) = 1024*x^11 - 2816x^9  + 2816x^7 - 1232x^5 + 220x^3 - 11x
    rts = find_zeros(T11, -1.0, 1.0)
    @test length(rts) == 11

    U9(x) = 512x^9 - 1024x^7 + 672x^5 - 160x^3 + 10x
    rts = find_zeros(U9, -1.0, 1.0)
    @test length(rts) == 9 # all rts real and in (-1,1) for Chebyshev

end


#http://www.chebfun.org/examples/roots/Tiger.html
tiger_tail(x) = 2*exp(.5*x) *(sin(5*x) + sin(101*x))
f1(x) = tiger_tail(x) - round(tiger_tail(x)) # (-2,1) then filter

f2(x) = (x-0.5)*(x-0.5001)*(x-4)*(x-4.05)*(x-9.3) # (0,10)
f3(x) = (x-3)^2*(x-4)^2 # (0,5)
delta1 = .01
f4(x) = (x-0.5)^2 * (x - (0.5+delta1))^2 # (0,1)
delta = .001
f5(x) = (x - 0.5)^3 * (x - (0.5+delta)) * (x - 1)
# f6(x) = (x - 0.5)^3 * (x - (0.5+delta))^3 * (x-4)*(x-(4+delta))*(x-4.2)^2

## test number of function calls with this
mutable struct CallableFunction
f
n
end
(F::CallableFunction)(x) = (F.n+=1; F.f(x))
f7 = CallableFunction(x -> sin(x-0.1), 0)


@testset "find_zeros harder ones" begin

    rts = find_zeros(f1, -2.0, 1.0, xatol=1e-4,xrtol=1e-4) # o/w misses a few
    rts = filter(u -> -1/4 < f1(u) < 1/4, rts)
    @test length(rts) == 345

    rts = find_zeros(f2, 0.0, 10.0)
    @test length(rts) == 5

    rts = find_zeros(f3, 0.0, 5.0)
    @test length(rts) == 2

    rts = find_zeros(f4, 0.0, 1.0)  # fussy
    @test length(rts) == 2

    # need to crank up xatol, xrtol to find smaller gaps (d=0.001)
    rts = find_zeros(f5, 0.0, 1.1, xatol=1e-5, xrtol=1e-5)
    @test length(rts) == 3

    # rts = find_zeros(f6, 0.0, 5.0, xatol=1e-5, xrtol=1e-5)
    # @test length(rts) == 5

    rts = find_zeros(f7, -pi/2, pi/2)
    @test f7.n <= 6000 # 3002

end
