# Root finding functions for Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMath.github.io/Roots.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMath.github.io/Roots.jl/dev)
[![Build Status](https://github.com/JuliaMath/Roots.jl/workflows/CI/badge.svg)](https://github.com/JuliaMath/Roots.jl/actions)
[![codecov](https://codecov.io/gh/JuliaMath/Roots.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/Roots.jl)

This package contains simple routines for finding roots of continuous
scalar functions of a single real variable. The `find_zero` function provides the
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
julia> using Roots

julia> f(x) = exp(x) - x^4;

julia> α₀,α₁,α₂ = -0.8155534188089607, 1.4296118247255556, 8.6131694564414;

julia> find_zero(f, (8,9), Bisection()) ≈ α₂ # a bisection method has the bracket specified
true

julia> find_zero(f, (-10, 0)) ≈ α₀ # Bisection is default if x in `find_zero(f,x)` is not a number
true


julia> find_zero(f, (-10, 0), Roots.A42()) ≈ α₀ # fewer function evaluations
true
```

For non-bracketing methods, the initial position is passed in as a
scalar:

```julia
julia> find_zero(f, 3) ≈ α₁  # find_zero(f, x0::Number) will use Order0()
true

julia> find_zero(f, 3, Order1()) ≈ α₁ # same answer, different method (secant)
true

julia> find_zero(sin, BigFloat(3.0), Order16()) ≈ π
true
```

The `find_zero` function can be used with callable objects:

```julia
julia> using Polynomials;

julia> x = variable();

julia> find_zero(x^5 - x - 1, 1.0) ≈ 1.1673039782614187
true
```

The function should respect the units of the `Unitful` package:

```julia
julia> using Unitful

julia> s, m  = u"s", u"m";

julia> g, v₀, y₀ = 9.8*m/s^2, 10m/s, 16m;


julia> y(t) = -g*t^2 + v₀*t + y₀
y (generic function with 1 method)

julia> find_zero(y, 1s)  ≈ 1.886053370668014s
true
```

Newton's method can be used without taking derivatives by hand. The
following use the `ForwardDiff` package:

```julia
julia> using ForwardDiff

julia> D(f) = x -> ForwardDiff.derivative(f,float(x))
D (generic function with 1 method)
```

Now we have:

```julia
julia> f(x) = x^3 - 2x - 5
f (generic function with 1 method)

julia> x0 = 2
2

julia> find_zero((f,D(f)), x0, Roots.Newton()) ≈ 2.0945514815423265
true
```

Automatic derivatives allow for easy solutions to finding critical
points of a function.

```julia
julia> using Statistics: mean, median

julia> as = rand(5);

julia> M(x) = sum([(x-a)^2 for a in as])
M (generic function with 1 method)

julia> find_zero(D(M), .5) ≈ mean(as)
true

julia> med(x) = sum([abs(x-a) for a in as])
med (generic function with 1 method)

julia> find_zero(D(med), (0, 1)) ≈ median(as)
true
```

### The CommonSolve interface

The
[DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
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

julia> solve(fx) ≈ 0.7728829591492101
true
```

With no default, and a single initial point specified, the default
`Order1` method is used.  The `solve` method allows other root-solving
methods to be passed, along with other options. For example, to use
the `Order2` method using a convergence criteria (see below) that
`|xₙ - xₙ₋₁| ≤ δ`, we could make this call:

```julia
julia> solve(fx, Order2(), atol=0.0, rtol=0.0) ≈ 0.7728829591492101
true
```

Unlike `find_zero`, which errors on non-convergence, `solve` returns
`NaN` on non-convergence.

This next example has a zero at `0.0`, but
for most initial values will escape towards `±∞`, sometimes causing a
relative tolerance to return a misleading value. Here we can see the
differences:

```julia
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

```julia
julia> f(x, p=2) = cos(x) - x/p
f (generic function with 2 methods)

julia> Z = ZeroProblem(f, pi/4)
ZeroProblem{typeof(f), Float64}(f, 0.7853981633974483)

julia> solve(Z, Order1()) ≈ 1.0298665293222586     # use p=2 default
true

julia> solve(Z, Order1(), p=3) ≈ 1.170120950002626 # use p=3
true
```

### Multiple zeros

The `find_zeros` function can be used to search for all zeros in a
specified interval. The basic algorithm essentially splits the interval into many
subintervals. For each, if there is a bracket, a bracketing algorithm
is used to identify a zero, otherwise a derivative free method is used
to search for zeros. This algorithm can miss zeros for various reasons, so the
results should be confirmed by other means.

```julia
julia> f(x) = exp(x) - x^4
f (generic function with 2 methods)

julia> find_zeros(f, -10,10) ≈ [α₀,α₁,α₂] # from above
true
```

The interval can also be specified using a structure with `extrema` defined, where `extrema` return two different values:

```julia
julia> using IntervalSets

julia> find_zeros(f, -10..10) ≈ [α₀,α₁,α₂]
true
```

(For tougher problems, the [IntervalRootFinding](https://github.com/JuliaIntervals/IntervalRootFinding.jl) package gives guaranteed results, rather than the heuristically identified values returned by `find_zeros`.)

### Convergence

For most algorithms, convergence is decided when

* The value `|f(x_n)| <= tol` with `tol = max(atol, abs(x_n)*rtol)`, or

* the values `x_n ≈ x_{n-1}` with tolerances `xatol` and `xrtol` *and*
  `f(x_n) ≈ 0` with a *relaxed* tolerance based on `atol` and `rtol`.

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

```julia
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

```julia
julia> f(x) = exp(x) - x^4
f (generic function with 2 methods)

julia> fzero(f, 8, 9) ≈ α₂   # bracketing
true

julia> fzero(f, -10, 0) ≈ α₀
true

julia> fzeros(f, -10, 10) ≈ [α₀, α₁, α₂]
true

julia> fzero(f, 3) ≈ α₁      # default is Order0()
true

julia> fzero(sin, big(3), order=16)  ≈ π # uses higher order method
true
```
