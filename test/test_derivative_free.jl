## tests of derivative free methods (orders 1,2,5,8,16)
using Roots
using Base.Test

orders = [1,2,5,8,16]

fn, x0, alpha = x -> cos(x) -x , 1.0, 0.7390851332151607
for order in orders
    @test_approx_eq fzero(fn, x0, order=order) alpha
    @test_approx_eq fzero(fn, big(x0), order=order) alpha
end


## failures
fn = x -> x^5 - x + 1 # initial steffensen step is a problem, so not secant, but others
x0 = -1.0
alpha = -1.1673039782614187

@test_approx_eq fzero(fn, x0, order=1) alpha
for order in orders[2:end]
    @test_throws Roots.ConvergenceFailed fzero(fn, x0, order=order)
end

fn = x -> x^20 - 1
x0, alpha = 0.5, 1.0
for order in orders[2:end]
    @test_throws Roots.ConvergenceFailed fzero(fn, x0, order=order)
end


## Issue near 0 that does not cross
fn = x -> cos(x) - 1
x0, alpha = .1, 0.0
#for order in orders[end]
#    @test_throws Roots.ConvergenceFailed fzero(fn, x0, order=order)
#end
for order in orders[end]
    ## test that these do not throw
    fn(fzero(fn, x0, order=order, ftol=1e-10))
end

## trivial
for order in orders
    @test_approx_eq fzero(x -> 0.0, 1, order=order)  1.0
    @test_approx_eq fzero(x -> x, 0.0, order=order)  0.0
end


## using hybrid mode: secant mode to get close, then 
fn = x -> x^5 - x + 1
x0, alpha = -1, -1.1673039782614187
@test_throws Roots.ConvergenceFailed fzero(fn, x0, order=8)
@test_approx_eq fzero(fn, fzero(fn, x0, order=1, ftol=0.25), order=8) alpha # and fewer steps order 1=25, hybrid=7 + 8
