# Simple tests to make sure Roots works with other
# packages. This is not part of runtests, as no expectation
# that CI should run these

using Roots
using Compat.Test

@testset "Test composability with other packages" begin

    orders = [Order0(), Order1(), Order2(), Order5(),
              Order8(), Order16()]
    
    # Unitful
    using Unitful
    @testset "find zero(s) with Unitful" begin
        
        
        s = u"s"; m = u"m"
        g = 9.8*m/s^2
        v0 = 10m/s
        y0 = 16m
        y(t) = -g*t^2 + v0*t + y0
        
        
        for order in orders
            @test find_zero(y, 1.8s, order) ≈  1.886053370668014s
        end
        
        @test find_zero(y, (1.8s, 1.9s), Bisection()) ≈ 1.886053370668014s
        
        xrts = find_zeros(y, 0s, 10s)
        @test length(xrts) == 1
        @test xrts[1] ≈  1.886053370668014s
    end
    
    
    # SymEngine
    using SymEngine
    @testset "find zero(s) with SymEngine" begin
        
        m = s = 1.0
        g = 9.8*m/s^2
        v0 = 10m/s
        y0 = 16m
        y(t) = -g*t^2 + v0*t + y0
        
        
        for order in orders
            @test find_zero(y, 1.8s, order) ≈  1.886053370668014s
        end
        
        @test find_zero(y, (1.8s, 1.9s), Bisection()) ≈ 1.886053370668014s
        
        xrts = find_zeros(y, 0s, 10s)
        @test length(xrts) == 1
        @test xrts[1] ≈  1.886053370668014s
    end



end
