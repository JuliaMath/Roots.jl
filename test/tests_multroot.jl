using Roots
using Polynomial
using Base.Test

x = Poly([1.0, 0.0])

p = (x-1)*(x-2)^2*(x-3)^3
zs, mults = fzero(p)
@test mults == [1,2,3]

p = (x-1)^2*(x-2)^2*(x-3)^4
zs, mults = fzero(p)
@test mults == [2,2,4]

p = (x-1)^2
zs, mults = multroot(p^14)
@test mults == [28]