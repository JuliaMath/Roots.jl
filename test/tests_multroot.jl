using Roots
using Polynomials
using Base.Test

## can use functions
f(x) = (x-1)*(x-2)^2*(x-3)^3
zs, mults = multroot(f)
@test mults == [1,2,3]

x = poly([0.0])

p = (x-1)*(x-2)^2*(x-3)^3
zs, mults = multroot(p)
@test mults == [1,2,3]

p = (x-1)^2*(x-2)^2*(x-3)^4
zs, mults = multroot(p)
@test mults == [2,2,4]

p = (x-1)^2
zs, mults = multroot(p^14)
@test mults == [28]


## test for real roots of polynomial functions

## for polynomials in Z[x], Q[x] can algorithm to be accurate for higher degree

@test fzeros(x -> x - 1)[1] == 1.0 # linear
f(x) = (x-1)*(x-2)*(x-3)^3*(x^2+1)
rts = fzeros(f)
@test maximum(abs(sort(rts) - [1.0, 2.0, 3.0])) <= 1e-12
x = poly([big(0)])
p = prod([x - i for i in 1:20])
Roots.real_roots(p) ## can find this
f(x) = (x-20)^5 - (x-20) + 1
a = fzeros(f)[1]
@assert abs(f(a)) <= 1e-14

## factor
factor(x -> (x-2)^4*(x-3)^9)
factor(x -> (x-1)^3 * (x-2)^3 * (x^5 - x + 1))
factor(x -> x*(x-1)*(x-2)*(x^2 + x + 1))

factor(x -> (x-1)^2 * (x-.99)^2 * (x-1.01)^2) ## can have issue with nearby roots (or high powers)

factor(x -> (x-1//1)^2 * (x-99//100)^2 * (x-101//100)^2) ## conversion is to Float, not Rational{Int}
factor(convert(Poly{Rational{Int}}, x -> (x-1//1)^2 * (x-99//100)^2 * (x-101//100)^2))






