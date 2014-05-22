using Roots
using Polynomial
using Base.Test

## can use functions
f(x) = (x-1)*(x-2)^2*(x-3)^3
zs, mults = multroot(f)
@test mults == [1,2,3]

x = Poly([1.0, 0.0])

p = (x-1)*(x-2)^2*(x-3)^3
zs, mults = multroot(p)
@test mults == [1,2,3]

p = (x-1)^2*(x-2)^2*(x-3)^4
zs, mults = multroot(p)
@test mults == [2,2,4]

p = (x-1)^2
zs, mults = multroot(p^14)
@test mults == [28]


## test for real roots
f(x) = (x-1)*(x-2)*(x-3)^3*(x^2+1)
rts = fzeros(f)
@test length(rts) == 3
