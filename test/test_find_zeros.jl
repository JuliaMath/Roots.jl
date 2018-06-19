using Roots
using Compat.Test



## test number of function calls with this
mutable struct CallableFunction
f
n
end
(F::CallableFunction)(x) = (F.n+=1; F.f(x))


#http://www.chebfun.org/examples/roots/Tiger.html
f1(x) = 2*exp(.5*x) *(sin(5*x) + sin(101*x))
tiger_tail(x) = f1(x) - round(f1(x)) # (-2,1) then filter

f2(x) = (x-0.5)*(x-0.5001)*(x-4)*(x-4.05)*(x-9.3) # (0,10)
f3(x) = (x-3)^2*(x-4)^2 # (0,5)
f4 = CallableFunction(x -> sin(x-0.1), 0)
delta = .001
f5(x) = (x - 0.5)^3 * (x - (0.5+delta)) * (x - 1)
f6(x) = (x - 0.5)^3 * (x - (0.5+delta))^3 * (x-4)*(x-(4+delta))*(x-4.2)^2

@testset "find_zeros" begin

    rts = find_zeros(tiger_tail, -2.0, 1.0)
    rts = filter(u -> -1/4 < tiger_tail(u) < 1/4, rts)
    @test length(rts) == 345

    rts = find_zeros(f2, 0.0, 10.0)
    @test length(rts) == 5

    rts = find_zeros(f3, 0.0, 5.0)
    @test length(rts) == 2

    rts = find_zeros(f4, -pi/2, pi/2)
    @test f4.n <= 10000  # 10,000 is *too* many

    rts = find_zeros(f5, 0.0, 1.5)
    @test length(rts) == 3

#    rts = find_zeros(f6, 0.0, 5.0)
#    @test length(rts) == 5
end
