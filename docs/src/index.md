# Roots.jl

Documentation for [Roots.jl](https://github.com/JuliaMath/Roots.jl)


## About

`Roots` is  a `Julia` package  for finding zeros of continuous
scalar functions of a single real variable using floating point numbers. That  is solving ``f(x)=0`` for ``x`` adjusting for floating-point idiosyncrasies.

The `find_zero` function provides the
primary interface. It supports various algorithms through the
specification of a method. These include:

* Bisection-like algorithms. For functions where a bracketing interval
  is known (one where ``f(a)`` and ``f(b)`` have alternate signs),
  there are several bracketing methods, including `Bisection`.  For
  most floating point number types, bisection occurs in a manner
  exploiting floating point storage conventions leading to an exact
  zero or a bracketing interval as small as floating point
  computations allows. Other methods include `A42`,
  `AlefeldPotraShi`, `Roots.Brent`, `Roots.Chandrapatlu`,
  `Roots.ITP`, `Roots.Ridders`, and ``12``-flavors of
  `FalsePosition`. The default bracketing method for
  the basic floating-point types is `Bisection` , as it is more robust to some inputs,
  but `A42` and `AlefeldPotraShi` typically converge in a few
  iterations and are more performant.


* Several derivative-free methods are implemented. These are specified
  through the methods `Order0`, `Order1` (the secant method), `Order2`
  (the Steffensen method), `Order5`, `Order8`, and `Order16`. The
  number indicates, roughly, the order of convergence. The `Order0`
  method is the default, and the most robust, as it finishes off with
  a bracketing method when a bracket is encountered, The higher order
  methods promise higher order (faster) convergence, though don't
  always yield results with fewer function calls than `Order1` or
  `Order2`. The methods `Roots.Order1B` and `Roots.Order2B` are
  superlinear and quadratically converging methods independent of the
  multiplicity of the zero.


* There are methods that require a derivative or two: `Roots.Newton`,
  `Roots.Halley` are classical ones, `Roots.QuadraticInverse`,
  `Roots.ChebyshevLike`, `Roots.SuperHalley` are others.
  `Roots.Schroder` provides a quadratic method, like Newton's method,
  which is independent of the multiplicity of the zero. The
  `Roots.ThukralXB`, `X=2`, `3`, `4`, or `5` are also multiplicity
  free. The `X` denotes the number of derivatives that need
  specifying. The `Roots.LithBoonkkampIJzerman{S,D}` methods remember
  `S` steps and use `D` derivatives.



## Basic usage

Consider  the polynomial   function  ``f(x) = x^5 - x + 1/2``. As a polynomial,  its roots, or  zeros, could  be identified with the  `roots` function of  the `Polynomials` package. However, even  that function uses a numeric method to identify   the values, as no  solution with radicals is available. That is, even for polynomials, non-linear root finders are needed to solve ``f(x)=0``. (Though polynomial root-finders can exploit certain properties not available for general non-linear functions.)

The `Roots` package provides a variety of algorithms for this  task. In this quick overview, only the  default ones  are illustrated.

For  the function ``f(x) = x^5 - x + 1/2`` a simple plot will show a zero  somewhere between ``-1.2`` and ``-1.0`` and two zeros near ``0.6``.

For the zero between two values at which the function changes sign, a
bracketing method is useful, as bracketing methods are guaranteed to
converge for continuous functions by the intermediate value
theorem. A bracketing algorithm will be used when the initial data is
passed as a tuple:

```jldoctest find_zero
julia> using Roots

julia> f(x) =  x^5 - x + 1/2
f (generic function with 1 method)

julia> find_zero(f, (-1.2,  -1)) ≈ -1.0983313019186336
true
```

The default algorithm is guaranteed to have an  answer nearly as accurate as is  possible  given the limitations of floating point  computations.

For the zeros "near" a point,  a non-bracketing method is often used, as generally  the algorithms are more efficient and can be  used in cases where a zero does  not cross the ``x`` axis. Passing just  the initial point will dispatch to  such a method:

```jldoctest find_zero
julia> find_zero(f,  0.6) ≈ 0.550606579334135
true
```


This finds  the answer  to the left of the starting point. To get the other nearby zero, a starting point closer to the answer can be used.

However,  an initial graph might convince one  that any of the up-to-``5`` real  roots  will   occur between ``-5``  and ``5``.  The `find_zeros` function uses  heuristics and a few of the  algorithms to  identify  all zeros between the specified range. Here  the method successfully identifies all  ``3``:

```jldoctest find_zero
julia> find_zeros(f, -5,  5)
3-element Vector{Float64}:
 -1.0983313019186334
  0.550606579334135
  0.7690997031778959
```
