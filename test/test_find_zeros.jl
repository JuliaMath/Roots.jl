using Roots
using Compat.Test


## test number of function calls with this
mutable struct CallableFunction
f
n
end
(F::CallableFunction)(x) = (F.n+=1; F.f(x))


@testset "find_zeros" begin

    xrts = find_zeros(x -> exp(x) - x^4, -5, 20)
    @test length(xrts) == 3

    n = 10
    xrts = find_zeros(x -> cos(x)/n, 0, 10)
    @test length(xrts) == 3

    W(n) = x -> prod(x-i for i in 1:n)
    xrts = find_zeros(W(20), -1, 21)
    @test length(xrts) == 20
 
        T11(x) = 1024*x^11 - 2816x^9  + 2816x^7 - 1232x^5 + 220x^3 - 11x
    xrts = find_zeros(T11, -1.0, 1.0)
    @test length(xrts) == 11 
    
    U9(x) = 512x^9 - 1024x^7 + 672x^5 - 160x^3 + 10x
    xrts = find_zeros(U9, -1.0, 1.0)
    @test length(xrts) == 9 # all xrts real and in (-1,1) for Chebyshev

    
    
    f1 = CallableFunction(x -> sin(x-0.1), 0)
    xrts = find_zeros(f1, -pi/2, pi/2)
    @test f1.n <= 6000 # ~ 3000
    
end



@testset "find_zeros harder ones" begin

    # non simple roots
    xrts = find_zeros(x -> (x-1)^2, 0, 2)
    @test length(xrts) == 1

    xrts = find_zeros(x->cos(x) + cos(2x), 0, 2pi)
    @test length(xrts) == 3


    
    #http://www.chebfun.org/examples/roots/Tiger.html
    tiger_tail(x) = 2*exp(.5*x) *(sin(5*x) + sin(101*x))
    f1(x) = tiger_tail(x) - round(tiger_tail(x)) # (-2,1) then filter
    
    f2(x) = (x-0.5)*(x-0.5001)*(x-4)*(x-4.05)*(x-9.3) # (0,10)
    f3(x) = (x-3)^2*(x-4)^2 # (0,5)
    delta1 = .01
    f4(x) = (x-0.5)^2 * (x - (0.5+delta1))^2 # (0,1)

    delta = .001
    f5(x) = (x - 0.5)^3 * (x - (0.5+delta)) * (x - 1)
    f6(x) = (x - 0.5)^3 * (x - (0.5+delta))^3 * (x-4)*(x-(4+delta))*(x-4.2)^2

    # These are fussy and depend on tolerances. We don't want
    # to have these run through CI, as it may fail on some machines
    # It would be nice were there an @test_allow_failure feature

    fussy_tests = Any[]
    
    xrts = find_zeros(f1, -2.0, 1.0, xatol=1e-4,xrtol=1e-4) # o/w misses a few
    xrts = filter(u -> -1/4 < f1(u) < 1/4, xrts)
    push!(fussy_tests, (:tiger_tail,length(xrts) == 345))

    xrts = find_zeros(f2, 0.0, 10.0, xatol=0.5*1e-4)
    push!(fussy_tests, (:f2, length(xrts) == 5))
    
    xrts = find_zeros(f3, 0.0, 5.0)
    push!(fussy_tests, (:f3, length(xrts) == 2))
    
    xrts = find_zeros(f4, 0.0, 1.0,xatol=0.005)  # delta1 < xatol
    push!(fussy_tests, (:f4, length(xrts) == 2))

    # need to crank up xatol, xrtol to find smaller gaps (d=0.001)
    xrts = find_zeros(f5, 0.0, 1.1, xatol=delta/2) 
    push!(fussy_tests, (:f5, length(xrts) == 3))

    xrts = find_zeros(f6, 0.0, 5.0, xatol=delta/2) # much harder due to 6th order near 0.5
    push!(fussy_tests, (:f6, length(xrts) == 5))


    W1(n) = x -> prod((x-i)^i for i in 1:n)
    xrts = find_zeros(W1(5), -1, 6) # o/w even powers hard
    push!(fussy_tests, (:W1, length(xrts) == 5))


    
    verbose=true
    if verbose
        println("")
        println("Testing fussy tests for find_zero")
        println("Test Summary:")
        print("Failed: "); [!b ? print(a) : print("") for (a,b) in fussy_tests]
        println("")
    end
    
end

