[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahub.com/docs/Roots/)
Linux: [![Build Status](https://travis-ci.org/JuliaMath/Roots.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Roots.jl)
Windows: [![Build status](https://ci.appveyor.com/api/projects/status/goteuptn5kypafyl?svg=true)](https://ci.appveyor.com/project/jverzani/roots-jl)

# Root finding functions for Julia

This package contains simple routines for finding roots of continuous
scalar functions of a single real variable. The `find_zero`function provides the
primary interface. It supports various algorithms through the
specification of a method. These include:

* Bisection-like algorithms. For functions where a bracketing interval
  is known (one where `f(a)` and `f(b)` have alternate signs), the
  `Bisection` method can be specified. For most floating point number
  types, bisection occurs in a manner exploiting floating point
  storage conventions. For others, an algorithm of Alefeld, Potra, and
  Shi is used. These methods are guaranteed to converge.

* Several derivative-free methods are implemented. These are specified
  through the methods `Order0`, `Order1` (the secant method), `Order2`
  (the Steffensen method), `Order5`, `Order8`, and `Order16`. The
  number indicates, roughly, the order of convergence. The `Order0`
  method is the default, and the most robust, but may take many more
  function calls to converge. The higher order methods promise higher
  order (faster) convergence, though don't always yield results with
  fewer function calls than `Order1` or `Order2`. The methods
  `Roots.Order1B` and `Roots.Order2B` are superlinear and quadratically converging
  methods independent of the multiplicity of the zero.

* There are historic methods that require a derivative or two:
  `Roots.Newton` and `Roots.Halley`.  `Roots.Schroder` provides a
  quadratic method, like Newton's method, which is independent of the
  multiplicity of the zero.

* There are several non-exported methods, such as, `Roots.Brent()`,
  `FalsePosition`, `Roots.A42`, `Roots.AlefeldPotraShi`,
  `Roots.LithBoonkkampIJzermanBracket`, and
  `Roots.LithBoonkkampIJzerman`.

Each method's documentation has additional detail.

Some examples:

```julia
using Roots
f(x) = exp(x) - x^4

# a bisection method has the bracket specified with a tuple, vector, or iterable
julia> find_zero(f, (8,9), Bisection())
8.613169456441398

julia> find_zero(f, (-10, 0))  # Bisection is default if x is a tuple and no method
-0.8155534188089606


julia> find_zero(f, (-10, 0), FalsePosition())  # just 11 function evaluations
-0.8155534188089607
```

For non-bracketing methods, the initial position is passed in as a
scalar:

```julia
## find_zero(f, x0::Number) will use Order0()
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
find_zero(x^5 - x - 1, 1.0)  # 1.1673039782614185
```

Or,

```julia
using Polynomials
x = variable()
find_zero(x^5 - x - 1, 1.0)  # 1.1673039782614185
```

The function should respect the units of the `Unitful` package:

```julia
using Unitful
s = u"s"; m = u"m"
g = 9.8*m/s^2
v0 = 10m/s
y0 = 16m
y(t) = -g*t^2 + v0*t + y0
find_zero(y, 1s)      # 1.886053370668014 s
```

Newton's method can be used without taking derivatives by hand. For example, if the
`ForwardDiff` package is available:

```julia
using ForwardDiff
D(f) = x -> ForwardDiff.derivative(f,float(x))
```

Now we have:

```
f(x) = x^3 - 2x - 5
x0 = 2
find_zero((f,D(f)), x0, Roots.Newton())   # 2.0945514815423265
```

Automatic derivatives allow for easy solutions to finding critical
points of a function.

```julia
## mean
using Statistics
as = rand(5)

function M(x)
  sum([(x-a)^2 for a in as])
end

find_zero(D(M), .5) - mean(as)	  # 0.0

## median
function m(x)
  sum([abs(x-a) for a in as])
end

find_zero(D(m), (0, 1)) - median(as)	# 0.0
```

### The CommonSolve interface

The
(DifferentialEquations)[https://github.com/SciML/DifferentialEquations.jl]
interface of setting up a problem; initializing the problem; then
solving the problem is also implemented using the methods
`ZeroProblem`, `init`, `solve!` and `solve`.

For example, we can solve a problem with many different methods, as follows:

```julia
julia> f(x) = exp(-x) - x^3
f (generic function with 1 method)

julia> x0 = 2.0
2.0

julia> fx = ZeroProblem(f,x0)
ZeroProblem{typeof(f), Float64}(f, 2.0)

julia> solve(fx)
0.7728829591492101
```

With no default, and a single initial point specified, the default
`Order1` method is used.  The `solve` method allows other root-solving
methods to be passed, along with other options. For example, to use
the `Order2` method using a convergence criteria (see below) that
`|xₙ - xₙ₋₁| ≤ δ`, we could make this call:

```julia
julia> solve(fx, Order2(), atol=0.0, rtol=0.0)
0.7728829591492101
```

Unlike `find_zero`, which errors on non-convergence, `solve` returns
`NaN` on non-convergence. This next example has a zero at `0.0`, but
for most initial values will escape towards `±∞`, sometimes causing a
relative tolerance to return a misleading value. Here we can see the
differences:

```julia
julia> f(x) = cbrt(x)*exp(-x^2)
f (generic function with 1 method)

julia> x0 = 0.1147
0.1147

julia> find_zero(f, 1.0, Roots.Order1()) # small relative value of f(xₙ), but not a mathematical zero
5.593607755898642

julia> find_zero(f, 1.0, Roots.Order1(), atol=0.0, rtol=0.0) # error as no check on `|f(xn)|`
ERROR: Roots.ConvergenceFailed("Stopped at: xn = 5.593607755898642. Too many steps taken. ")
[...]

julia> fx = ZeroProblem(f, x0);

julia> solve(fx, Roots.Order1(), atol=0.0, rtol=0.0) # NaN, not an error
NaN

julia> fx = ZeroProblem((f, D(f)), x0) # higher order methods can identify zero of this function
ZeroProblem{Tuple{typeof(f), var"#1#2"{typeof(f)}}, Float64}((f, var"#1#2"{typeof(f)}(f)), 0.1147)

julia> solve(fx, Roots.LithBoonkkampIJzerman(2,1), atol=0.0, rtol=0.0)
0.0
```

The iterator interface can be useful for hybrid zero-finding algorithms. Here we illustrate the iterator produced by `init` on a rapidly convergent example:

```
julia> fx = ZeroProblem(sin, (3,4))
ZeroProblem{typeof(sin), Tuple{Int64, Int64}}(sin, (3, 4))

julia> p = init(fx, Roots.A42())
Roots.A42: x₀: [3.0, 3.5]

julia> 

julia> for _ ∈ p; @show p; end
p = x₁: [3.14156188905231, 3.1416247553172956]

p = x₂: [3.141592653589793, 3.1415926535897936]

julia> p
Converged to x₂: 3.141592653589793
```

### Multiple zeros

The `find_zeros` function can be used to search for all zeros in a
specified interval. The basic algorithm essentially splits the interval into many
subintervals. For each, if there is a bracket, a bracketing algorithm
is used to identify a zero, otherwise a derivative free method is used
to search for zeros. This algorithm can miss zeros for various reasons, so the
results should be confirmed by other means.

```julia
f(x) = exp(x) - x^4
find_zeros(f, -10, 10)  # -0.815553…,  1.42961…,  8.61317…
```

(For tougher problems, the [IntervalRootFinding](https://github.com/JuliaIntervals/IntervalRootFinding.jl) package gives guaranteed results, rather than the heuristically identified values returned by `find_zeros`.)

### Convergence

For most algorithms, convergence is decided when

* The value |f(x_n)| < tol with `tol = max(atol, abs(x_n)*rtol)`, or

* the values x_n ≈ x_{n-1} with tolerances `xatol` and `xrtol` *and*
  f(x_n) ≈ 0 with a *relaxed* tolerance based on `atol` and `rtol`.

The algorithm stops if

* it encounters an `NaN` or an `Inf`, or

* the number of iterations exceed `maxevals`, or

* the number of function calls exceeds `maxfnevals`.

If the algorithm stops and the relaxed convergence criteria is met,
the suspected zero is returned. Otherwise an error is thrown
indicating no convergence. To adjust the tolerances, `find_zero`
accepts keyword arguments `atol`, `rtol`, `xatol`, and `xrtol`, as
seen in some examples above.


The `Bisection` and `Roots.A42` methods are guaranteed to converge
even if the tolerances are set to zero, so these are the
defaults. Non-zero values for `xatol` and `xrtol` can be specified to
reduce the number of function calls when lower precision is required.

```
julia> fx = ZeroProblem(sin, (3,4));

julia> p = init(fx, Bisection(), xatol=1/4)
Bisection: x₀: [3.0, 4.0]

julia> for _ ∈ p; @show p; end
p = x₁: [3.0, 3.5]
p = x₂: [3.0, 3.25]
p = x₃: [3.125, 3.25]
```

## An alternate interface

This functionality is provided by the `fzero` function, familiar to
MATLAB users. `Roots` also provides this alternative interface:

* `fzero(f, x0::Real; order=0)` calls a
  derivative-free method. with the order specifying one of
  `Order0`, `Order1`, etc.

* `fzero(f, a::Real, b::Real)` calls the `find_zero` algorithm with the
  `Bisection` method.

* `fzeros(f, a::Real, b::Real)` will call `find_zeros`.

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
fzero(sin, big(3), order=16)  # 3.141592653589793...
```

### Technical difference between find_zero and fzero

The `fzero` function is not identical to `find_zero`. When a function, `f`,
is passed to `find_zero` the code is specialized to the function `f`
which means the first use of `f` will be slower due to compilation,
but subsequent uses will be faster. For `fzero`, the code is not
specialized to the function `f`, so the story is reversed.
