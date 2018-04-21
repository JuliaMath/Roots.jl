[![Roots](http://pkg.julialang.org/badges/Roots_0.6.svg)](http://pkg.julialang.org/?pkg=Roots&ver=0.6)  
Linux: [![Build Status](https://travis-ci.org/JuliaMath/Roots.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Roots.jl)
Windows: [![Build status](https://ci.appveyor.com/api/projects/status/goteuptn5kypafyl?svg=true)](https://ci.appveyor.com/project/jverzani/roots-jl)

# Root finding functions for Julia


This package contains simple routines for finding roots of continuous
scalar functions of a single real variable. The `find_zero`function provides the
primary interface. It support various algorithms through the
specification of an algorithm. These include:

* Bisection-like algorithms. For functions where a bracketing interval
  is known (one where $f(a)$ and $f(b)$ have alternate signs), the
  `Bisection` method can be specified with a guaranteed
  convergence. For most floating point number types, bisection occurs
  in a manner exploiting floating point storage conventions. For
  others, an algorithm of Alefeld, Potra, and Shi is used.

  For typically faster convergence -- though not guaranteed -- the
  `FalsePosition` method can be specified. This methods has one of 12
  implementations which implement a modified secant method to
  accelerate convergence.

* Several derivative-free methods are implemented. These are specified
  through the methods `Order0`,
  `Order1` (the secant method), `Order2` (the Steffensen method),
  `Order5`, `Order8`, and `Order16`. The number indicates roughly the
  order of convergence. The `Order0` method is  default, and the most robust, but
  may take many more function calls. The higher order ones promise
  higer order convergence, though don't always yield results with fewer
  function calls.

* There are two historic methods that require a derivative:
  `Roots.Newton` and `Roots.Halley`. (Neither is currently exported.)
  If a derivative is not given, an automatic derivative is found using
  the `ForwardDiff` package.

Some examples: 


```
# a bisection method has the bracket specified with a tuple or vector
julia> find_zero(f, (8,9), Bisection())
8.613169456441398

julia> find_zero(f, (-10, 0))  # Bisection if x is a tuple and no method
-0.8155534188089606


julia> find_zero(f, (-10, 0), FalsePosition())  # just 19 function evaluations
-0.8155534188089607
	
julia> find_zero(f, 3)         # default is Order0()
1.4296118247255556

julia> find_zero(f, 3, Order1()) # same answer, different method
1.4296118247255556

julia> find_zero(sin, BigFloat(3.0), Order16())
3.141592653589793238462643383279502884197169399375105820974944592307816406286198
```


The `find_zero` function can be used with callable objects:

```julia
using SymEngine
@vars x
find_zero(x^5 - x - 1, 1.0)
```

Or,

```julia
using Polynomials
x = variable(Int)
fzero(x^5 - x - 1, 1.0)
```

The function should respect the units of the `Unitful` package:

```
using Unitful
s = u"s"; m = u"m"
g = 9.8*m/s^2
v0 = 10m/s
y0 = 16m
y(t) = -g*t^2 + v0*t + y0
find_zero(y, 1s)
```

Newton's method can be used without taking derivatives:

```
f(x) = x^3 - 2x - 5
x0 = 2
find_zero(f, x0, Roots.Newton())
```

Automatic derivatives allow for easy solutions to finding critical
points of a function.

```julia
## mean
as = rand(5)
function M(x) 
  sum([(x-a)^2 for a in as])
end
fzero(D(M), .5) - mean(as)	  # 0.0

## median
function m(x) 
  sum([abs(x-a) for a in as])

end
fzero(D(m), 0, 1)  - median(as)	# 0.0
```

### Multiple zeros

The `find_zeros` function can be used to search for all zeros in a
specified interval. The basic algorithm splits the interval into many
subintervals. For each, if there is a bracket a bracketing algorithm
is used to identify a zero, otherwise a derivative free method is used
to check. This algorithm can miss zeros for various reasons, so the
results should be confirmed by other means.

```
f(x) = exp(x) - x^4
find_zeros(f, -10, 10)
```


### Convergence

For most algorithms (besides the `Bisection` ones) convergence is decided when

* The value f(x_n) ≈ 0 with tolerances `atol` and `rtol` *or*

* the values x_n ≈ x_{n-1} with tolerances `xatol` and `xrtol` *and*
f(x_n) ≈ 0 with a *relaxed* tolerance based on `atol` and `rtol`.

* an algorithm encounters an `NaN` or `Inf` and yet f(x_n) ≈ 0 with a *relaxed* tolerance based on `atol` and `rtol`.

There is no convergence if the number of iterations exceed `maxevals`,
or the number of function calls exceeds `maxfnevals`.

The tolerances may need to be adjusted. To determine if convergence
occurs due to f(x_n) ≈ 0, it is necessary to consider that even if
`xstar` is the correct answer mathematically, due to floating point
roundoff it is expected that f(xstar) ≈ f'(xstar) ⋅ xstar ⋅ ϵ. A
relative error accounts for the value of `x`, but the default
tolerance may need adjustment if the derivative is large near the
zero, as the default is a bit aggressive. On the other hand, the
absolute tolerance might seem too relaxed. 

To determine if convergence is determined as x_n ≈ x_{n-1} the check
on f(x_n) ≈ 0 is done as algorithms can be fooled by asymptotes, or
other areas where the tangent lines have large slopes.

The `Bisection` and `Roots.A42` methods will converge, so the tolerances are ignored.

## An alternate interface

For MATLAB users, this functionality is provided by the `fzero`
function. `Roots` also provides this alternative.

The function `fzero` which dispatches to an appropriate
algorithm based on its argument(s):

* `fzero(f, a::Real, b::Real)` and `fzero(f,
  bracket::Vector)` call the `find_zero` algorithm with the
  `Bisection` method.
  
* `fzero(f, x0::Real; order::Int=0)` calls a
  derivative-free method. with the order specified matching one of
  `Order0`, `Order1`, etc.
  
* `fzeros(f, a::Real, b::Real; no_pts::Int=200)` will call `find_zeros`.

* The function `secant_method`, `newton`, and `halley` provide direct
  access to those methods.


## Usage examples

```julia
f(x) = exp(x) - x^4
## bracketing
fzero(f, 8, 9)		          # 8.613169456441398
fzero(f, -10, 0)		      # -0.8155534188089606
fzeros(f, -10, 10)            # -0.815553, 1.42961  and 8.61317 

## use a derivative free method
fzero(f, 3)			          # 1.4296118247255558

## use a different order
fzero(sin, 3, order=16)		  # 3.141592653589793
```





----

Some additional documentation can be read [here](http://nbviewer.ipython.org/url/github.com/JuliaLang/Roots.jl/blob/master/doc/roots.ipynb?create=1).
