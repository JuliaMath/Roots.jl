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
only requires one new  function call per  step  so  can be very effective. Often  function evaluations are the  slowest part of  the computation and, as  well, no derivative is  needed. Because  it  can be  very efficient, the secant  method  is used in  the default method  of `find_zero` when  called with a single initial starting point.

[Steffensen's](https://en.wikipedia.org/wiki/Steffensen%27s_method) method is a quadratically converging. derivative-free method  which uses a secant  line  based on ``x_n`` and ``x_n + f(x_n)``.  Though of  higher  order, it requires  additional function calls per step and depends on a  good initial starting value. Other  derivative free methods are available, trading off  increased function calls for higher-order convergence. They may be  of interest when arbitrary  precision is needed. A  measure of efficiency is ``q^{1/r}`` where ``q`` is the order of convergence and ``r`` the number of function calls per step.   With this measure, the secant method  would be ``\approx (1.618)^{1/1}`` and Steffensen's  would be less (``2^{1/2}``).

----

```@docs
Secant
Steffensen
Order5
Order8
Order16
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
Roots.Esser
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

| Type            | Method                       | Order                  | F evals | Asymptotic efficiency                 |
|:--------------- | :--------------------------- | :--------------------- | :------ | :------------------------------------ |
| Hybrid          | Order0                       |                        |         | ``\approx 1.618\dots``                |
| Derivative Free | Secant                       | ``\varphi=1.618\dots`` | ``1``   | ``1.618\dots``                        |
| Derivative Free | Steffensen                   | ``2``                  | ``2``   | ``1.414\dots``                        |
| Derivative Free | Order5                       | ``5``                  | ``4``   | ``1.495\dots``                        |
| Derivative Free | Order8                       | ``8``                  | ``4``   | ``1.681\dots``                        |
| Derivative Free | Order16                      | ``16``                 | ``5``   | ``1.718\dots``                        |
| Classical       | Newton                       | ``2``                  | ``2``   | ``1.414\dots``                        |
| Classical       | Halley                       | ``3``                  | ``3``   | ``1.442\dots``                        |
| Classical       | QuadraticInverse             | ``3``                  | ``3``   | ``1.442\dots``                        |
| Classical       | ChebyshevLike                | ``3``                  | ``3``   | ``1.442\dots``                        |
| Classical       | SuperHalley                  | ``3``                  | ``3``   | ``1.442\dots``                        |
| MultiStep       | LithBoonkkampIJzerman{S,D}   | ``p^s=\sum p^k(d+\sigma_k)`` | ``D+1`` | varies, ``1.92\dots`` max       |
| Bracketing      | BisectionExact               | ``1``                  | ``1``   | ``1``                                 |
| Bracketing      | A42                          | ``(2 + 7^{1/2})``      | ``3,4`` |``(2 + 7^{1/2})^{1/3} = 1.6686\dots``  |
| Bracketing      | AlefeldPotraShi              |                        | ``3,4`` | ``1.618\dots``                        |
| Bracketing      | Brent                        | ``\leq 1.89\dots``     | ``1``   | ``\leq 1.89\dots``                    |
| Bracketing      | ITP                          | ``\leq \varphi``         | ``1``   | ``\leq \varphi``                      |
| Bracketing      | Ridders                      | ``1.83\dots``          | ``2``   | ``1.225\dots``                          |
| Bracketing      | FalsePosition                | ``1.442\dots``         | ``1``   | ``1.442\dots``                        |
| Bracketing      | LithBoonkkampIJzermanBracket | ``2.91``               | ``3``   | ``1.427\dots``                        |
| Robust          | King                         | ``\varphi=1.618\dots`` | ``2``   | ``1.272\dots``                        |
| Robust          | Esser                        | ``2``                  | ``3``   | ``1.259\dots``                        |
| Robust          | Schroder                     | ``2``                  | ``3``   | ``1.259\dots``                        |
| Robust          | Thukral3                     | ``3``                  | ``4``   | ``1.316\dots``                        |
| Robust          | Thukral4                     | ``4``                  | ``5``   | ``1.319\dots``                        |
| Robust          | Thukral5                     | ``5``                  | ``6``   | ``1.307\dots``                        |



## Convergence

Identifying when an algorithm converges or diverges requires specifications of tolerances  and convergence criteria.

In the case of exact bisection, convergence is mathematically
guaranteed. For floating point numbers, either an *exact* zero is
found, or the bracketing interval can be subdivided into ``[a_n,b_n]``
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
not near a ``0``. As such, for non-bracketing methods, a check on the
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

The initial naming scheme used `fzero` instead  of `fzeros`, following the name of the  MATLAB function [fzero](https://www.mathworks.com/help/matlab/ref/fzero.html). This interface  is not recommended, but, for now, still maintained.

```@docs
Roots.fzero
```

## Tracking iterations

It is possible to add the keyword argument `verbose=true`  when calling the `find_zero` function to get detailed information about the solution and data from each iteration. To save this data a `Tracks` object may be passed in to `tracks`.

----

```@docs
Roots.Tracks
```
