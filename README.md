# Root finding functions for Julia

This package contains simple routines for finding roots of continuous
scalar functions of a single real variable. The basic interface is
through the function `fzero` which dispatches to an appropriate
algorithm based on its arguments:

* `fzero(p::Union(Function, Poly), a::Real, b::Real)` and
  `fzero(p::Union(Function, Poly), bracket::Vector)` call the
  `find_zero` algorithm to find a root within the bracket `[a,b]`. For convenience,
  `fzeros(f::Function, a::Real, b::Real)` will split the interval `[a,b]` and apply `fzero` to each
  bracketing subintervals to search for all zeros in the intervals

* `fzero(p::Union(Function, Poly), x0::Real; order::Int=0)` calls a
  derivative-free method. The default is slow but more robust to the
  quality of the initial guess. If a bracket can be found, the
  returned value will satisfy `f(x)==0.0` or
  `f(prevfloat(x))*f(nextfloat(x)) <= 0`. Otherwise a value is returned
  with `f(x)` very nearly 0.0 or an error is thrown. For faster
  convergence and less memory usage, an order can be
  specified. Possible values are 1, 2, 5, 8, and 16. The order 2
  Steffensen method can be the fastest, but is more sensitive to a good
  initial guess. The order 8 method is more robust and often as
  fast. The higher-order method may be faster when using `Big` values.

* `fzero(p::Union(Function, Poly), x0::Real, bracket::Vector)` calls
  a derivative-free algorithm with initial guess `x0` with steps constrained
  to remain in the specified bracket.

* `fzero(p::Poly)` calls `multroot`, An improvement on the `roots`
  function of the `Polynomial` package when multiple roots are
  present. Follows algorithms due to Zeng, ["Computing multiple roots
  of inexact polynomials", Math. Comp. 74 (2005),
  869-903](http://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01692-8/home.html).


For historical purposes, there are implementations of Newton's method
(`newton`), Halley's method (`halley`), and the secant method
(`secant_method`). For the first two, if derivatives are not
specified, they will be computed using the `PowerSeries` package.


## Usage examples

```
f(x) = exp(x) - x^4
## bracketing
fzero(f, [8, 9])		# 8.613169456441398
fzero(f, -10, 0)		# -0.8155534188089606
fzeros(f, -10, 10)		# -0.815553, 1.42961  and 8.61317 

## use a derivative free method
fzero(f, 3)			# 1.4296118247255558

## use a different order
fzero(sin, 3, order=16)		# 3.141592653589793

## BigFloat values yield more precision
fzero(sin, BigFloat(3.0))	# 3.1415926535897932384...with 256 bits of precision
```

Polynomial root finding is a bit better when multiple roots are present.

```
using Roots, Polynomial
x = Poly([1,0])
p = (x-1)^2 * (x-3)
fzero(p)			# compare to roots(p)
fzero(p, [2.5, 3.5])
```

The well-known methods can be used with or without supplied
derivatives. If not specified, the `PowerSeries` package is used for
automatic differentiation.

```
f(x) = exp(x) - x^4
fp(x) = exp(x) - 4x^3
fpp(x) = exp(x) - 12x^2
newton(f, fp, 8)		# 8.613169456441398
newton(f, 8)	
halley(f, fp, fpp, 8)
halley(f, 8)
secant_method(f, 8, 8.5)
```

The automatic derivatives allow for easy solutions to finding critical
points of a function.

```
## mean
as = rand(5)
function M(x) 
  sum([(x-a)^2 for a in as])
end
fzero(D(M), .5) - mean(as)	# 0.0

## median
function m(x) 
  sum([abs(x-a) for a in as])
end
fzero(D(m), 0, 1)  - median(as)	# 0.0
```
