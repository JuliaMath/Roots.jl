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

```jldoctest readme
julia> using Roots

julia> f(x) = exp(x) - x^4;

julia> find_zero(f, (8,9), Bisection()) # a bisection method has the bracket specified
8.6131694564414

julia> find_zero(f, (-10, 0))  # Bisection is default if x in `find_zero(f,x)` is not a number
-0.8155534188089607


julia> find_zero(f, (-10, 0), Roots.A42())  # fewer function evaluations
-0.8155534188089607
```

For non-bracketing methods, the initial position is passed in as a
scalar:

```jldoctest readme
julia> find_zero(f, 3)   # find_zero(f, x0::Number) will use Order0()
1.4296118247255556

julia> find_zero(f, 3, Order1()) # same answer, different method (secant)
1.4296118247255556

julia> find_zero(sin, BigFloat(3.0), Order16())
3.141592653589793238462643383279502884197169399375105820974944592307816406286198
```

The `find_zero` function can be used with callable objects:

```jldoctest readme
julia> import Pkg; Pkg.add("Polynomials"); # if not installed
[...]

julia> using Polynomials;
[...]

julia> x = variable();

julia> find_zero(x^5 - x - 1, 1.0)
1.1673039782614187
```

The function should respect the units of the `Unitful` package:

```jldoctest readme
julia> import Pkg; Pkg.add("Unitful"); # if not installed
[...]

julia> using Unitful

julia> s = u"s"; m = u"m"
m

julia> g = 9.8*m/s^2
9.8 m s⁻²

julia> v0 = 10m/s
10 m s⁻¹

julia> y0 = 16m
16 m

julia> y(t) = -g*t^2 + v0*t + y0
y (generic function with 1 method)

julia> find_zero(y, 1s)      # 1.886053370668014 s
1.8860533706680143 s
```

Newton's method can be used without taking derivatives by hand. The
following use the `ForwardDiff` package:

```jldoctest readme
julia> import Pkg; Pkg.add("ForwardDiff"); # if not installed
[...]

julia> using ForwardDiff

julia> D(f) = x -> ForwardDiff.derivative(f,float(x))
D (generic function with 1 method)
```

Now we have:

```jldoctest readme
julia> f(x) = x^3 - 2x - 5
f (generic function with 1 method)

julia> x0 = 2
2

julia> find_zero((f,D(f)), x0, Roots.Newton())
2.0945514815423265
```

Automatic derivatives allow for easy solutions to finding critical
points of a function.

```jldoctest readme
julia> using Statistics: mean, median

julia> as = rand(5);

julia> M(x) = sum([(x-a)^2 for a in as])
M (generic function with 1 method)

julia> find_zero(D(M), .5) - mean(as)
0.0

julia> med(x) = sum([abs(x-a) for a in as])
med (generic function with 1 method)

julia> find_zero(D(med), (0, 1)) - median(as)
0.0
```

### The CommonSolve interface

The
(DifferentialEquations)[https://github.com/SciML/DifferentialEquations.jl]
interface of setting up a problem; initializing the problem; then
solving the problem is also implemented using the methods
`ZeroProblem`, `init`, `solve!` and `solve`.

For example, we can solve a problem with many different methods, as follows:

```jldoctest readme
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

```jldoctest readme
julia> solve(fx, Order2(), atol=0.0, rtol=0.0)
0.7728829591492101
```

Unlike `find_zero`, which errors on non-convergence, `solve` returns
`NaN` on non-convergence. This next example has a zero at `0.0`, but
for most initial values will escape towards `±∞`, sometimes causing a
relative tolerance to return a misleading value. Here we can see the
differences:

```jldoctest readme
julia> f(x) = cbrt(x)*exp(-x^2)
f (generic function with 1 method)

julia> x0 = 0.1147
0.1147

julia> find_zero(f, 1.0, Roots.Order1()) # stopped as |f(xₙ)| ≤ |xₙ|ϵ
5.53043654482315

julia> find_zero(f, 1.0, Roots.Order1(), atol=0.0, rtol=0.0) # error as no check on `|f(xn)|`
ERROR: Roots.ConvergenceFailed("Algorithm failed to converge")
[...]

julia> fx = ZeroProblem(f, x0);

julia> solve(fx, Roots.Order1(), atol=0.0, rtol=0.0) # NaN, not an error
NaN

julia> fx = ZeroProblem((f, D(f)), x0); # higher order methods can identify zero of this function

julia> solve(fx, Roots.LithBoonkkampIJzerman(2,1), atol=0.0, rtol=0.0)
0.0
```

Functions may be parameterized, as illustrated:

```jldoctest readme
julia> f(x, p=2) = cos(x) - x/p
f (generic function with 2 methods)

julia> Z = ZeroProblem(f, pi/4)
ZeroProblem{typeof(f), Float64}(f, 0.7853981633974483)

julia> solve(Z, Order1()) # use p=2 default
1.0298665293222586

julia> solve(Z, Order1(), p=3)
1.170120950002626
```

### Multiple zeros

The `find_zeros` function can be used to search for all zeros in a
specified interval. The basic algorithm essentially splits the interval into many
subintervals. For each, if there is a bracket, a bracketing algorithm
is used to identify a zero, otherwise a derivative free method is used
to search for zeros. This algorithm can miss zeros for various reasons, so the
results should be confirmed by other means.

```jldoctest readme
julia> f(x) = exp(x) - x^4
f (generic function with 2 methods)

julia> find_zeros(f, -10, 10)  # -0.815553…,  1.42961…,  8.61317…
3-element Vector{Float64}:
 -0.8155534188089606
  1.4296118247255556
  8.613169456441398
```

The interval can also be specified using a structure with `extrema` defined, where `extrema` return two different values:

```
julia> import Pkg; Pkg.add("IntervalSets");
[...]

julia> using IntervalSets
[...]

julia> find_zeros(f, -10..10)
3-element Vector{Float64}:
 -0.8155534188089606
  1.4296118247255556
  8.613169456441398
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

If the algorithm stops and the relaxed convergence criteria is met,
the suspected zero is returned. Otherwise an error is thrown
indicating no convergence. To adjust the tolerances, `find_zero`
accepts keyword arguments `atol`, `rtol`, `xatol`, and `xrtol`, as
seen in some examples above.


The `Bisection` and `Roots.A42` methods are guaranteed to converge
even if the tolerances are set to zero, so these are the
defaults. Non-zero values for `xatol` and `xrtol` can be specified to
reduce the number of function calls when lower precision is required.

```jldoctest readme
julia> fx = ZeroProblem(sin, (3,4));

julia> solve(fx, Bisection(); xatol=1/16)
3.125
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

### Usage examples

```jldoctest readme
julia> f(x) = exp(x) - x^4
f (generic function with 2 methods)

julia> fzero(f, 8, 9)    # bracketing
8.6131694564414

julia> fzero(f, -10, 0)
-0.8155534188089607

julia> fzeros(f, -10, 10)
3-element Vector{Float64}:
 -0.8155534188089606
  1.4296118247255556
  8.613169456441398

julia> fzero(f, 3)       # default is Order1()
1.4296118247255556

julia> fzero(sin, big(3), order=16)  # uses higher order method
3.141592653589793238462643383279502884197169399375105820974944592307816406286198
```
