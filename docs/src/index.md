# Roots.jl

Documentation for [Roots.jl](https://github.com/JuliaMath/Roots.jl)


## About

`Roots` is  a `Julia` package  for finding zeros of continuous
scalar functions of a single real variable using floating point numbers. That  is solving ``f(x)=0`` for ``x`` adjusting for floating-point idiosyncrasies.

The `find_zero` function provides the
primary interface. It supports various algorithms through the
specification of a method. These include:

* Bisection-like methods. For functions where a bracketing interval
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


* Several derivative-free methods. These are specified
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


* Methods requiring one or more derivatives: `Roots.Newton`,
  `Roots.Halley` are classical ones, `Roots.QuadraticInverse`,
  `Roots.ChebyshevLike`, `Roots.SuperHalley` are others.
  `Roots.Schroder` provides a quadratic method, like Newton's method,
  which is independent of the multiplicity of the zero. The
  `Roots.ThukralXB`, `X=2`, `3`, `4`, or `5` are also multiplicity
  free. The `X` denotes the number of derivatives that need
  specifying. The `Roots.LithBoonkkampIJzerman{S,D}` methods remember
  `S` steps and use `D` derivatives.
