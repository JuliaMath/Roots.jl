# Reference/API

The `Roots`  package provides several  different  algorithms  to solve `f(x)=0`.


```@index
Pages = ["reference.md"]
```

```@setup reference
using Roots
```

```@meta
DocTestSetup = quote
  using Roots
end
```

```@meta
CurrentModule = Roots
```

```@docs
Roots
```

##  The `find_zero`  and `find_zeros` functions

There are  two main  functions:  `find_zero`   to  identify  a  zero  of  ``f``  given  some initial starting  value or  bracketing interval and  `find_zeros` to heuristically identify  all  zeros in  a specified interval.


```@docs
find_zero
find_zeros
```

## CommonSolve interface

The problem-algorithm-solve interface is a pattern popularized in `Julia` by the `DifferentialEquations.jl` suite of packages. This can be used as an alternative to `find_zero`. Unlike `find_zero`, `solve` will return `NaN` on non-convergence.

```@docs
Roots.solve!
Roots.solve
Roots.ZeroProblem
```



##  Classical  methods  based on derivatives

We begin  by  describing  the classical methods even though they are not necessarily  recommended  because they require more work of the  user,  as they give insight into  why there  are a variety  of methods available.

The classical  methods of [Newton](https://en.wikipedia.org/wiki/Newton%27s_method) and  [Halley](https://en.wikipedia.org/wiki/Halley%27s_method) utilize information about  the function  and  its derivative(s) in  an  iterative manner  to converge to  a zero of  ``f(x)`` given an initial starting value.

Newton's method is   easily described:

From  an initial point,  the  next  point  in  the iterative algorithm is found by identifying the  intersection of  the ``x``    axis  with  the tangent line of ``f`` at the initial  point. This is repeated until convergence  or the realization that   convergence won't happen for the  initial point. Mathematically,

``x_{n+1}  =  x_{n}  - f(x_n)/f'(x_n).``

Some facts  are helpful  to  understand the different methods  available in `Roots`:

* For Newton's method there is a formula for the error: Set
  ``\epsilon_n = \alpha - x_n``, where ``\alpha`` is the zero, then
  ``\epsilon_{n+1} = -f''(\xi_n)/(2f'(\xi_n) \cdot \epsilon_n^2,``
  here ``\xi_n`` is some value between ``\alpha`` and ``x_n``.

* The error term, when of the form ``|\epsilon_{n+1}| \leq
  C\cdot|\epsilon_n|^2``, can be used to identify an interval around
  ``\alpha`` for which convergence is guaranteed. Such convergence is
  termed *quadratic* (order 2).  For floating point solutions,
  quadratic convergence and a well chosen initial point can lead to
  convergence in 4 or 5 iterations. In general, convergence is termed
  order ``q`` when ``|\epsilon_{n+1}| \approx C\cdot|\epsilon_n|^q``

* The term ``-f''(\xi_n)/(2f'(\xi_n)`` indicates possible issues  when ``f''``  is  too big  near ``\alpha``  or  ``f'`` is too small  near ``\alpha``. In particular if ``f'(\alpha)  =  0``, there need  not be quadratic  convergence, and convergence  can   take many  iterations. A  zero   for which ``f(x) = (x-\alpha)^{1+\beta}\cdot g(x)``, with ``g(\alpha) \neq 0``  is called *simple* when ``\beta=0`` and  non-simple when  ``\beta >  0``. Newton's method is quadratic near *simple  zeros* and need not be quadratic  near  *non-simple* zeros.
As well,  if  ``f''`` is too  big near ``\alpha``, or  ``f'`` too small near  ``\alpha``, or ``x_n``  too  far  from  ``\alpha`` (that is,  ``|\epsilon_n|>1``) the  error  might actually increase and convergence is not guaranteed.

* The explicit form of  the error function can  be used to guarantee convergence for functions with a certain shape (monotonic, convex functions where the sign of ``f''`` and ``f'`` don't change). Quadratic convergence may only occur once the algorithm is near the zero.

* The number of function evaluations  per step for Newton's method is 2.

----

```@docs
Roots.Newton
Roots.Halley
Roots.QuadraticInverse
Roots.ChebyshevLike
Roots.SuperHalley
```

Newton and Halley's method are members of this family of methods:

```@docs
Roots.LithBoonkkampIJzerman{S,D}
```


## Derivative free methods

The [secant](https://en.wikipedia.org/wiki/Secant_method) method replaces the  derivative term in Newton's method with the slope of a secant line using two prior values:

``x_{n+1} = x_n - (\frac{f(x_n)-f(x_{n-1})}{x_n - x_{n-1}})^{-1}\cdot  f(x_n).``

Though the secant  method   has  convergence  rate of  order ``\approx 1.618`` -- i.e., is not quadratic --  it
only requires one new  function call per  step  so  can be very effective. Often  function evaluations are the  slowest part of  the computation and, as  well, no derivative is  needed. Because  it  can be  very efficient, the secant  method  is used in  the default method  of `find_zero` when  called with a single initial starting point. The `Roots.Sidi` methods generalize the secant method.

[Steffensen's](https://en.wikipedia.org/wiki/Steffensen%27s_method) method is a quadratically converging. derivative-free method  which uses a secant  line  based on ``x_n`` and ``x_n + f(x_n)``.  Though of  higher  order, it requires  additional function calls per step and depends on a  good initial starting value. Other  derivative free methods are available, trading off  increased function calls for higher-order convergence. They may be  of interest when arbitrary  precision is needed. A  measure of efficiency is ``q^{1/r}`` where ``q`` is the order of convergence and ``r`` the number of function calls per step.   With this measure, the secant method  would be ``\approx (1.618)^{1/1}`` and Steffensen's  would be less (``2^{1/2}``).

----

```@docs
Secant
Order1
Steffensen
Order2
Order5
Order8
Order16
Roots.Sidi
Roots.Muller
```


## Bracketing methods

The [bisection](https://en.wikipedia.org/wiki/Bisection_method) method identifies a zero of a *continuous* function between ``a`` and ``b``  when  ``f(a)`` and  ``f(b)`` have different  signs. (The interval ``[a,b]`` is called a bracketing interval when ``f(a)\cdot  f(b)  <0``.)  The basic  algorithm is particularly simple, an interval  ``[a_i,b_i]`` is  split  at  ``c =  (a_i+b_i)/2``. Either  ``f(c)=0``,  or one  of  ``[a_i,c]``  or  ``[c,b_i]`` is a bracketing  interval,  which is  called  ``[a_{i+1},b_{i+1}]``. From this  description,  we  see  that  ``[a_i,b_i]`` has length  ``2^{-i}`` times the length of ``[a_0,b_0]``, so  the intervals will eventually terminate by finding  a zero, ``c``,  or converge  to a zero. This convergence is slow (the efficiency  is only ``1``, but guaranteed. For  `16`-, `32`-, and `64`-bit  floating point  values, a  reinterpretation  of  how the midpoint  (``c``) is found  leads  to convergence  in  no more  than   ``64`` iterations, unlike the midpoint found above, where some cases can take many more steps to converge.

In floating point,  by  guaranteed  convergence we have either an exact zero or a bracketing interval  consisting   of  two  adjacent floating point values. When applied to *non*-continuous  functions,  this algorithm  will identify   an exact  zero or  a zero crossing   of the function. (E.g., applied  to  ``f(x)=1/x`` it  will  find  ``0``.)

The default selection of  midpoint described above includes no information  about the function ``f`` beyond its  sign. Algorithms exploiting  the shape of the function  can be significantly more efficient. For example, the bracketing method `Roots.AlefeldPotraShi` due to [Alefeld, Potra, and Shi](https://dl.acm.org/doi/10.1145/210089.210111) has  efficiency ``\approx 1.6686``. This method  is  also   used in the  default method for `find_zero` when a  single initial starting point is given if a bracketing interval is identified.

----

```@docs
Bisection
Roots.A42
Roots.AlefeldPotraShi
Roots.Brent
Roots.Chandrapatla
Roots.Ridders
Roots.ITP
Roots.ModAB
Roots.RegulaFalsi
FalsePosition
Roots.LithBoonkkampIJzermanBracket
Roots.BracketedHalley
Roots.BracketedChebyshev
Roots.BracketedSchroder
```

## Non-simple zeros

The order of convergence for most methods is for *simple* zeros, values ``\alpha`` where ``f(x) = (x-\alpha) \cdot g(x)``, with ``g(\alpha)`` being non-zero. For methods which are of order ``k`` for non-simple zeros, usually an additional function call is needed per step. For example, this is the case for `Roots.Newton` as compared to `Roots.Schroder`.

Derivative-free methods for non-simple zeros have the following implemented:

```@docs
Roots.King
Roots.Order1B
Roots.Esser
Roots.Order2B
```


For non-simple zeros, Schroder showed an additional derivative can  be used to yield quadratic convergence based on Newton's method:

```@docs
Roots.Schroder
```


A family of methods for non-simple zeros which require ``k`` derivatives to be order ``k``, with ``k=2`` yielding Schroder's method, are implemented in:

```@docs
Roots.AbstractThukralBMethod
```



## Hybrid  methods

A useful  strategy  is   to  begin with a non-bracketing  method and switch to a bracketing method should a bracket be encountered. This  allows   for the identification of zeros which are not surrounded by a bracket, and have guaranteed convergence  should a bracket be  encountered.  It  is  used  by default by `find_zero(f,a)`.

```@docs
Roots.Order0
```

## All zeros

The `find_zeros` function heuristically scans an interval for all zeros using a combination of bracketing and non-bracketing methods. The `AllZeros` method may be passed to `solve` to call this.

```@docs
Roots.AllZeros
```

## Rates of convergence

The order of a method is ``q``, where ``e_{i+1} \approx
e_i^q``. Newton's method is famously quadratic **for** simple roots;
the secant method of order ``q \approx \varphi=1.618\dots``. However,
``p=2`` calls are needed for Newton's method, and only ``p=1`` for the
secant method. The asymptotic efficiency is ``q^{1/p}``, which
penalizes function calls. There are other order ``k`` methods taking
``k`` function calls per step, e.g., Halley's; others take fewer, as
seen below. Many use inverse quadratic steps, others inverse
cubic--these have order ``q`` solving ``q^{s+1}-2q^s+1`` (``s=3`` for
quadratic). For robust methods, generally ``1`` additional function
call is needed to achieve the convergence rate, `Schroder` being a
good example.

| Type            | Method                       | Order                  | F evals | Asymptotic efficiency                 | Reference |
|:--------------- | :--------------------------- | :--------------------- | :------ | :------------------------------------ |:----------|
| Hybrid          | Order0                       |                        |         | ``\approx 1.618\dots``                |[Order0](@cite)|
| Derivative Free | Secant                       | ``\varphi=1.618\dots`` | ``1``   | ``1.618\dots``                        ||
| Derivative Free | Steffensen                   | ``2``                  | ``2``   | ``1.414\dots``                        ||
| Derivative Free | Order5                       | ``5``                  | ``4``   | ``1.495\dots``                        |[Order5](@cite)|
| Derivative Free | Order8                       | ``8``                  | ``4``   | ``1.681\dots``                        |[Order8](@cite)|
| Derivative Free | Order16                      | ``16``                 | ``5``   | ``1.718\dots``                        |[Order16](@cite)|
| Derivative Free | Muller                       | ``1.839\dots``         | ``1``   | ``1.839\dots``                        |[Muller's method](https://en.wikipedia.org/wiki/Muller%27s_method)|
| Classical       | Newton                       | ``2``                  | ``2``   | ``1.414\dots``                        ||
| Classical       | Halley                       | ``3``                  | ``3``   | ``1.442\dots``                        ||
| Classical       | QuadraticInverse             | ``3``                  | ``3``   | ``1.442\dots``                        |[Chebyshev](@cite)|
| Classical       | ChebyshevLike                | ``3``                  | ``3``   | ``1.442\dots``                        |[ChebyshevLike](@cite)|
| Classical       | SuperHalley                  | ``3``                  | ``3``   | ``1.442\dots``                        |[SuperHalley](@cite)|
| MultiStep       | LithBoonkkampIJzerman{S,D}   | ``p^s=\sum p^k(d+\sigma_k)`` | ``D+1`` | varies, ``1.92\dots`` max       |[LithBoonkkampIJzerman](@cite)|
| Bracketing      | BisectionExact               | ``1``                  | ``1``   | ``1``                                 ||
| Bracketing      | A42                          | ``(2 + 7^{1/2})``      | ``3,4`` |``(2 + 7^{1/2})^{1/3} = 1.6686\dots``  |[AlefeldPotraShi](@cite)|
| Bracketing      | AlefeldPotraShi              |                        | ``3,4`` | ``1.618\dots``                        |[AlefeldPotraShi](@cite)|
| Bracketing      | Brent                        | ``\leq 1.89\dots``     | ``1``   | ``\leq 1.89\dots``                    ||
| Bracketing      | ITP                          | ``\leq \varphi``       | ``1``   | ``\leq \varphi``                      |[ITP](@cite)|
| Bracketing      | Ridders                      | ``1.83\dots``          | ``2``   | ``1.225\dots``                          |[Ridders](@cite)|
| Bracketing      | RegulaFalsi{:classic}       | ``1``                  | ``1``   | ``1``                                 |[RegulaFalsi](@cite)|
| Bracketing      | RegulaFalsi{:Illinois}      | ``1.442\dots``         | ``1``   | ``1.442\dots``                        |[RegulaFalsi](@cite)|
| Bracketing      | RegulaFalsi{:AndersonBjork}  | ``1.681\dots``         | ``1``   | ``1.681\dots``                        |[RegulaFalsi](@cite)|
| Bracketing      | RegulaFalsi{:Ford4}          | ``1.681\dots``         | ``1``   | ``1.681\dots``                        |[RegulaFalsi](@cite)|
| Bracketing      | ModAB                        | ``≈1.7\dots``          | ``1``   | ``1.7\dots``                          |[ModAB](@cite)|
| Bracketing      | LithBoonkkampIJzermanBracket | ``2.91``               | ``3``   | ``1.427\dots``                        |[LithBoonkkampIJzerman](@cite)|
| Robust          | King                         | ``\varphi=1.618\dots`` | ``2``   | ``1.272\dots``                        |[King](@cite)|
| Robust          | Esser                        | ``2``                  | ``3``   | ``1.259\dots``                        |[Esser](@cite)|
| Robust          | Schroder                     | ``2``                  | ``3``   | ``1.259\dots``                        |[Schroder](@cite)|
| Robust          | Thukral3B                     | ``3``                  | ``4``   | ``1.316\dots``                        |[ThukralB](@cite)|
| Robust          | Thukral4B                     | ``4``                  | ``5``   | ``1.319\dots``                        |[ThukralB](@cite)|
| Robust          | Thukral5B                     | ``5``                  | ``6``   | ``1.307\dots``                        |[ThukralB](@cite)|



## Convergence

Identifying when an algorithm converges or diverges requires specifications of tolerances  and convergence criteria.

In the case of exact bisection, convergence is mathematically
guaranteed. For floating point numbers, either an *exact* zero is
found, or the bracketing interval can *not* be subdivided into ``[a_n,b_n]``
with ``a_n`` and ``b_n`` being adjacent floating point values. That is
``b_n-a_n`` is as small as possible in floating point numbers. This can
be considered a stopping criteria in ``\Delta x``. For early termination
(less precision but fewer function calls) a tolerance can be given so
that if ``\Delta_n=b_n-a_n`` is small enough the algorithm stops
successfully.  In floating point, assessing if ``b_n \approx a_n``
requires two tolerances: a *relative* tolerance, as the minimal
differences in floating point values depend on the size of ``b_n`` and
``a_n``, and an absolute tolerance for values near ``0``. The values
`xrtol` and `xatol` are passed to the `Base.isapprox` function to
determine closeness.

Relying on the closeness of two ``x`` values will not be adequate for
all problems, as there are examples where the difference
``\Delta_n=|x_n-x_{n-1}|`` can be quite small, ``0`` even, yet ``f(x_n)`` is
not near a ``0``. As well, a final step which should make $\Delta_n$ small enough, might get derailed due to floating point issues. As such, for non-bracketing methods and for some non-strict bracketing methods, a check on the
size of ``f(x_n)`` is also used. As we find floating point
approximations to ``\alpha``, the zero, we must consider values small
when ``f(\alpha(1+\epsilon))`` is small. By Taylor's approximation, we
can expect this to be around
``\alpha\cdot \epsilon \cdot f'(\alpha)``.
That is, small depends on the size of ``\alpha`` and the
derivative at ``\alpha``.  The former is handled by both relative and absolute
tolerances (`rtol` and `atol`).  The size of ``f'(\alpha)`` is problem
dependent, and can be accommodated by larger relative or absolute
tolerances.

When an algorithm returns  an  `NaN` value,  it terminates. This  can  happen near convergence or  may indicate some issues.  Early termination is checked for convergence  in the  size  of ``f(x_n)`` with a relaxed tolerance when `strict=false` is specified (the default).

!!! note "Relative tolerances  and assessing  `f(x) ≈ 0`"
    The use of  relative tolerances  to  check  if   ``f(x)  \approx  0`` can lead  to spurious  answers  where  ``x`` is very large   (and  hence the relative  tolerance  is large). The return of  very  large solutions  should  be checked against expectations  of the  answer.


Deciding if an algorithm won't  terminate is  done  through  counting the number or  iterations performed; the default  adjusted through `maxiters`. As  most  algorithms are superlinear, convergence happens rapidly near  the answer, but  all the algorithms  can take  a while  to  get near  an  answer, even   when progress  is made. As  such, the maximum must be large enough to consider linear cases, yet small enough to avoid too many steps when an algorithm is non-convergent.


Convergence criteria are method dependent and are determined by  the  `Roots.assess_convergence`  methods.


```@docs
Roots.assess_convergence
```

Default tolerances  are specified through the `Roots.default_tolerances` methods.

```@docs
Roots.default_tolerances
```




## Simplified versions


The abstractions and many checks for  convergence employed by `find_zero` have a performance cost. When that is a critical concern, there are  several "simple" methods provided which can offer improved performance.

```@docs
Roots.secant_method
Roots.bisection
Roots.muller
Roots.newton
Roots.dfree
```


## MATLAB interface

The initial naming scheme used `fzero` instead  of `find_zero`, following the name of the  MATLAB function [fzero](https://www.mathworks.com/help/matlab/ref/fzero.html). This interface  is not recommended, but, for now, still maintained.

```@docs
fzero
fzeros
```

## Tracking iterations

To get detailed information about the solution and data from each iteration  a `Tracks` object may be passed in to `tracks`.

----

```@docs
Roots.Tracks
```

## Implementation details

The basic shell to find zeros of a function allowing for different methods is a bit more complicated than just programming a simple loop for each method. `Julia`'s multiple dispatch makes the approach quite uniform, allowing `find_zero(f, x0, M)` to solve ``f(x) = 0`` in many different manners with the same iteration scheme, described below.

### Arguments to `find_zero`

The `find_zero` function is called with three necessary things:

* A function (or functions)
* A starting point (or points)
* A method (or methods); the basic unit of dispatch

#### `Callable_Function`

The function or functions needed for evaluation is typically just some callable object. However, for some methods, like Newton's method, two or more functions are needed. When two (or more) functions are needed they can be passed in a tuple (e.g. `(sin,cos)`) *or* as a single function returning `(f, f/f', f'/f'', ...)`. To keep a similar interface, the `Callable_Function` struct is given.

Consider these examples:

```@example details
using Roots # hide
F = Roots.Callable_Function(Roots.Secant(), sin)
F(3)
```

or

```@example details
F = Roots.Callable_Function(Roots.Newton(), (sin, cos))
F(3)
```

or

```@example details
F = Roots.Callable_Function(Roots.Newton(), x -> (sin(x), sin(x)/cos(x)))
F(3)
```

The details of how many values are needed for the algorithm is provided by `Roots.fn_argout(M)`:

```@example details
Roots.fn_argout(Roots.Secant()), Roots.fn_argout(Roots.Newton())
```

#### The initial guess

For Newton's method, the initial guess is a nearby single number; for the Secant method, the initial guess is a nearby set of two numbers; for Bisection, the initial guess is two numbers which form a *bracketing* interval.

Eventually, all ``x`` values specified are coerced to a floating point value.

For a single value, there isn't more to discuss; for two values, there is more detail.

For a secant step, if only a single number is specified, another one is generated. The utility function `_default_secant_step` does this. The suggested step depends on the scale of the single value.

```@example details
Roots._default_secant_step(1.0), Roots._default_secant_step(1e16)
```

For a bracketing method, `Roots.assert_bracket` is called on ``f(a)`` and ``f(b)``.

```@example details
Roots.assert_bracket(-1.0, 1.0) # checks if product of their signs is negative
```

Bracketing intervals can be specified by any object that has `extrema` defined and for which that call returns two distinct values.


#### Methods

Methods are specified through an instance, such as `Roots.Secant()`. Many are exported and do not need qualification. There are some methods which depend on parameters. Examples are

```@example details
Roots.LithBoonkkampIJzerman(3,2), Roots.Sidi(2)
```

As mentioned elsewhere, when a method is omitted there are two different defaults:

* if the initial guess has distinct values returned by `extrema`, these are assumed to be a bracketing interval and a bracketing algorithm is called (which depends on the number type, though typically this is bisection over 64-bit floating point values.)

* if the initial guess is a single number, then the `Order0` method is called. This is a *hybrid* method which begins with *guarded* secant steps but will switch to a bracketing method *if* a bracketing interval is identified. Bracketing methods are guaranteed to converge. Other hybrid methods can be specified by specifying *two* methods, e.g.,  `find_zero(f, x0, M, N)`.


### Setup

The main algorithm uses the `CommonSolve` setup where a problem is initialized and then solved.

#### ZeroProblem

The `init` function of `CommonSolve` must dispatch on an object. `ZeroProblem` bundles the function(s) and initial guess into such an object.

```@example details
Roots.ZeroProblem(sin, 4)
```

#### `init`


This `ZeroProblem` object and the method(s) are then passed to `init` which

* Uses the method and function(s) to create a `Callable_Function` object, possibly incorporating a parameter passed to `find_zero` or `solve`.

* Uses the method, the `Callable_Function` object, and the initial guess to create a "state" object.

##### State objects

The basic state object holds the last two ``x``-values and their corresponding ``f(x)`` values. Some methods require more detail and so there are specialized states. An example would be the `ModABState` which also records a count and a flag for an algorithm which switches during the iteration. The state object has ``x``-values of type `T <: AbstractFloat`; the ``f(x)`` values may be of a different type.

```@example details
M = Roots.Secant()
f, x0 = sin, 3.0
Z = ZeroProblem(f, x0)
F = Roots.Callable_Function(M, Z.F)
state = Roots.init_state(M, F, Z.x₀)
```

##### Options object

The method, the state, and the keyword arguments are used to create an options object which holds the tolerances. Keyword arguments override the default arguments which are method and state dependent (the state carries the types).

```@example details
o = Roots.init_options(M, state; atol=1e-8)
```

The options record the absolute and relative tolerance in the ``x`` and ``f(x)`` values (`atol` above sets `abstol`--the fourth argument--in the struct), the maximum number of iterations and two flags for convergence checking.



### Algorithm

With the method, the function, the state object, the options object, and a logging object, the algorithm is ready to be run.

There are two main functions:

#### `update_state`

This method depends on the method, the function, the state, the options, and the logger. The method updates the state according to the algorithm (States are immutable, so it really creates a new state each time). This method typically generates the next ``x`` value and its ``y`` value and then shifts these to the most recent. The `update_state` method passes back the updated state object *and* a flag in cases some numeric value is identified which should stop the algorithm.

#### `assess_convergence`

Convergence is assessed differently based on the method. Bracketing methods may only check on the size of ``x_n - x_{n-1}``, where as secant-like methods might also check the size of the residual ``f(x_n)``. If these are small, the function returns a flag to stop.


#### `decide_convergence`

The two function above are in a loop, the number of times each is called is compared and compared to the `maxiters` value for the method (or passed through `maxiters`). If too many steps are needed the algorithm stops.

A stopped algorithm might or might not have *converged*.

It is possible the value is numerically close to an answer, but the tolerances might have been too tight or the algorithm might just have failed. Deciding between these scenarios is done by `decide_convergence` which (if instructed by the value of `strict`) utilizes a relaxed set of tolerances to see if convergence has occurred. Most bracketing methods don't use this, as convergence upto floating point limits can typically be guaranteed, but secant-like methods utilize this.

By allowing very strict tolerance in the algorithm with relaxed tolerance in the final decision, more accuracy can be expected than with just a relaxed initial tolerance.

If convergence is decided, then the value identified is returned.

If convergence is not decided, then `solve` returns a typed `NaN` value. The `find_zero` method errors when this happens, `solve` just returns `NaN`.


## Bibliography

```@bibliography
```
