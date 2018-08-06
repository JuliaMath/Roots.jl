using Compat.Test
import Roots.fzero

## Test `fzero` interface to `find_zero`
## test `fzeros` interface for functions

orders = [0,1,2,5,8,16]

### derivative free
fn, xstar, x0, br = x -> sin(x) - x/2, 1.89549426703398094714, 2.0, [pi/2, pi]
@test fzero(fn, x0)  ≈ xstar
for o in orders
    @test fzero(fn, x0, order=o)  ≈ xstar
end
### bisection
@test fzero(fn, br)  ≈ xstar
@test fzero(x->5-x,5,10)==5

### test tolerances
fn, xstar, x0, br = x -> x^5 - x - 1, 1.1673039782614187, 1.0, [1.0, 2.0]
@test fzero(fn, x0, order=1)  ≈ xstar

@test_throws Roots.ConvergenceFailed fzero(fn, x0, order=1, maxevals=2)
#@test abs(fzero(fn, x0, order=1, ftol=1e-2) - xstar) > 1e-5
#@test abs(fzero(fn, x0, order=1, xtol=1e-2) - xstar) > 1e-10

@test abs(find_zero(fn, br, Roots.A42(), xatol=1e-1) - xstar) > 1e-5


## various tests
## issue #29, basically f(a) or f(b) so big we get NaN for guess from bracket.
f =  x -> x + exp(x)
@test fzero(f, [-1e6, 1e6])  ≈ -0.5671432904097838

f =  x -> 1/x - 1
@test fzero(f, [0, 2])  ≈ 1.0

## test infinite range
@test fzero(x -> x, [-Inf, Inf])  ≈ 0.0

##################################################
## fzeros function
rts = 1:5
@test all((abs.(fzeros(x -> prod([x-r for r in rts]),0,10)) .- collect(1:5)) .<= 1e-15)


fn = x -> cos(10*pi*x)
@test length(fzeros(fn, 0, 1)) == 10

### issue with fzeros and roots near 'b'
@test 0 <  maximum(fzeros(x -> sin(x) - 1/1000*x, 0, pi)) < pi

