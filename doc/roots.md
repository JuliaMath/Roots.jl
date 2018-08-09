# The Roots package

The `Roots` package contains simple routines for finding zeros of
continuous scalar functions of a single real variable.  A zero of $f$
is a value $c$ where $f(c) = 0$.  The basic interface is through the
function `find_zero`, which through multiple dispatch can handle many different cases.

In the following, we will use `Plots` for plotting and `ForwardDiff`
to take derivatives.

```
using Roots
using Plots, ForwardDiff; pyplot();
```

## Bracketing

For a function $f$ (univariate, real-valued) a *bracket* is a pair $ a < b $
for which $f(a) \cdot f(b) < 0$. That is the function values have
different signs at $a$ and $b$. If
$f$ is a continuous function this ensures
([Bolzano](https://en.wikipedia.org/wiki/Intermediate_value_theorem))
there will be a zero in the interval $[a,b]$.  If $f$ is not
continuous, then there must be a point $c$ in $[a,b]$ where the function
"jumps" over $0$.

Using bisection, such values can be found, up to floating point
roundoff. That is, given `f(a) * f(b) < 0`, a value `c` with `a < c < b` can be
found where either `f(c) == 0.0` or  `f(prevfloat(c)) * f(c) < 0` or
`f(c) * f(nextfloat(c)) < 0`.

To illustrate, consider the function $f(x) = \cos(x) - x$. From the
graph we see readily that $[0,1]$ is a bracket (which we emphasize
with an overlay):

```figure
f(x) = cos(x) - x
plot(f, -2, 2)
plot!([0,1], [0,0], linewidth=2)
```

The `Roots` package includes the bisection algorithm through
`find_zero`. We use a vector or tuple to specify the initial condition
and `Bisection()` to specify the algorithm:

```
x = find_zero(f, (0, 1), Bisection())   
x, f(x)
```

For this function we see that `f(x)` is `0.0`.

----

Next consider $f(x) = \sin(x)$. A known zero is $\pi$. Trigonometry
tells us that $[\pi/2, 3\pi/2]$ will be a bracket. In this call `Bisection()`
is not specified, as it will be the default when the initial value is
specified as a pair of numbers:

```
f(x) = sin(x)
x = find_zero(f, (pi/2, 3pi/2))
x, f(x)
```

This value of `x` does not exactly produce a zero, however, it is as close as can be:

``` f(prevfloat(x)) * f(x) < 0.0 || f(x) * f(nextfloat(x)) < 0.0 ```

That is, at `x` the function is changing sign.

From a mathematical perspective, a zero is guaranteed for a
*continuous* function. However, the computer algorithm doesn't assume
continuity, it just looks for changes of sign. As such, the algorithm
will  identify discontinuities, not just zeros. For example:

```
find_zero(x -> 1/x, (-1, 1))
```


The endpoints and function values can even be infinite:

```
find_zero(x -> Inf*sign(x), (-Inf, Inf))  # Float64 only
```


The basic algorithm used for bracketing when the values are simple
floating point values is a modification of the bisection method. For big float values, an
algorithm due to Alefeld, Potra, and Shi is used.

```
find_zero(sin, (big(3), big(4)))    # uses a different algorithm then for (3,4)
```

By default, bisection will converge to machine tolerance. This may
provide more accuracy than desired. A tolerance may be specified to
terminate early, thereby utilizing fewer resources. For example, this
uses 19 steps to reach accuracy to $10^{-6}$ (without specifying `xatol` it uses
51 steps):

```
rt = find_zero(sin, (3.0, 4.0), xatol=1e-6)
rt - pi
```

There are more efficient algorithms than `Bisection` implemented:
`Roots.A42()`, `Roots.AlefeldPotraShi()`, `Roots.Brent()`, and 12
flavors of `Roots.FalsePosition()`. The first few have guaranteed
convergence and use, generally, fewer function calls than bisection
(for sufficiently smooth functions). The false position methods are
also faster, though may fail to converge on certain problems.


## Non-bracketing problems

Bracketing methods may have guaranteed convergence, but in general may
require more function calls than are needed to produce an answer.  If
a good initial guess is known, then the `find_zero` function provides
an interface to some different iterative algorithms that are more
efficient for good initial guesses. Unlike bracketing methods, these
algorithms may not converge to the desired root if the initial guess
is not well chosen.

The default algorithm is modeled after an algorithm used for
[HP-34 calculators](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf). This
algorithm is designed to be more forgiving of the quality of the
initial guess at the cost of possibly performing many more steps than
other algorithms, as if the algorithm encounters a bracket, a
bracketing algorithm with guaranteed convergence will be used.

For example, the answer to our initial problem is visibly seen from a
graph to be near 1. Given this,
the zero is found through:

```
f(x) = cos(x) - x
x = find_zero(f , 1)
x, f(x)
```

For the polynomial $f(x) = x^3 - 2x - 5$, an initial guess of 2 seems reasonable:

```
f(x) = x^3 - 2x - 5
x = find_zero(f, 2)
x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
```

For even more precision, `BigFloat` numbers can be used

```
x = find_zero(sin, big(3))
x, sin(x), x - pi
```

### Higher order methods

The default call to `find_zero` uses the first order secant method and then
possibly a bracketing method, which involves potentially many more function
calls. There may be times where a more efficient algorithm is
sought.  For such, a higher-order method might be better suited. There
are algorithms `Order1` (secant method), `Order2`
([Steffensen](http://en.wikipedia.org/wiki/Steffensen's_method)),
`Order5`, `Order8`, and `Order16`. The order 1 or 2 methods are
generally quite efficient. The even higher order ones are potentially
useful when more precision is used. These algorithms are accessed by
specifying the method after the initial starting point:

```
f(x) = 2x - exp(-x)
x = find_zero(f, 1, Order1())      
x, f(x)
```


```
f(x) = (x + 3) * (x - 1)^2
x = find_zero(f, -2, Order2())
x, f(x)
```


```
x = find_zero(f, 2, Order8())
x, f(x)
```

The latter shows that zeros need not be simple zeros  to be found.
A simple zero, $c$,has $f(x) = (x-c) \cdot g(x)$ where $g(c) \neq 0$.
Generally speaking, non-simple zeros are
expected to take many more function calls, as the methods are no
longer super-linear. This is the case here, where `Order2` uses $55$
function calls, `Order8` uses $41$, and `Order0` takes, a comparable,
$42$, as it does not identify a bracketing interval.

To investigate an algorithm and its convergence, the argument
`verbose=true` may be specified.


For some functions, adjusting the default tolerances may be necessary
to achieve convergence. The tolerances include `atol` and `rtol`, which are
used to check if $f(x_n) \approx 0$;
`xatol` and `xrtol`, to check if $x_n \approx x_{n-1}$; and `maxevals` and `maxfnevals` to limit the
number of steps in the algorithm or function calls.



## Classical methods

The package provides some classical methods for root finding:
`Roots.newton`, `Roots.halley`, and `Roots.secant_method`. (Currently
these are not exported, so must be prefixed with the package name to
be used.) We can see how each works on a problem studied by Newton
himself. Newton's method uses the function and its derivative:

```
f(x) = x^3 - 2x - 5
fp(x) = 3x^2 - 2
x = Roots.newton(f, fp, 2)
x, f(x)
```

To see the algorithm in progress, the argument `verbose=true` may be
specified.

Alternatively, `Roots.Newton()` can be specified as the method for `find_zero`. The
functions are specified using a tuple:

```
find_zero((f,fp), 2, Roots.Newton())
```

The secant method typically needs two starting points, though a second
one is computed if only one is give. Here we start with 2 and 3,
specified through a tuple:

```
x = Roots.secant_method(f, (2,3))
x, f(x)
```

Starting with a single point is also supported:

```
Roots.secant_method(f, 2)
```

(This is like `Order1()`, but the implementation is significantly
faster, as the framework is bypassed, and fewer checks on convergence
are used. This method can be used when speed is very important.)

Halley's method has cubic convergence, as compared to Newton's
quadratic convergence. It uses the second derivative as well:

```
fpp(x) = 6x
x = Roots.halley(f, fp, fpp, 2)
x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
```

(Halley's method takes 3 steps, Newton's 4, but Newton's uses 5
function calls to Halley's 10.)

For many functions, their derivatives can be computed automatically. The
`ForwardDiff` package provides a means. Here we define an operator `D`
to compute a derivative:

```
using ForwardDiff
D(f) = x -> ForwardDiff.derivative(f, float(x))
D(f, n) = n > 1 ? D(D(f),n-1) : D(f)
```


```
Roots.newton(f, D(f), 2)    
```

Or, for Halley's method:

```
Roots.halley(f, D(f), D(f,2), 2)  
```


## Finding critical points

The `D` function, defined above, makes it straightforward to find critical points
(typically where the derivative is $0$ but also where it is undefined). For example, the critical
point of the function $f(x) = 1/x^2 + x^3, x > 0$ near $1.0$ is where
the derivative is $0$ and can be found through:

```
f(x) = 1/x^2 + x^3
find_zero(D(f), 1)
```

For more complicated expressions, `D` will not work, and other means
of finding a derivative can be employed. In
this example, we have a function that models the flight
of an arrow on a windy day:

```
function flight(x, theta)
 	 k = 1/2
	 a = 200*cosd(theta)
	 b = 32/k
	 tand(theta)*x + (b/a)*x - b*log(a/(a-x))
end
```



The total distance flown is when `flight(x) == 0.0` for some `x > 0`:
This can be solved for different `theta` with `find_zero`. In the
following, we note that `log(a/(a-x))` will have an asymptote at `a`,
so we start our search at `a-5`:

```
function howfar(theta)
	 a = 200*cosd(theta)
	 find_zero(x -> flight(x, theta), a-5)
end
```	 

To visualize the trajectory if shot at 45 degrees, we have:


```figure
theta = 45
plot(x -> flight(x,  theta), 0, howfar(theta))
```


To maximize the range we solve for the lone critical point of `howfar`
within reasonable starting points. In this example, the derivative can not be taken automatically with
`D`. So,  here we use a central-difference approximation and start the
search at 45 degrees--the angle which maximizes the trajectory on a
non-windy day:

```
h = 1e-5
howfarp(theta) = (howfar(theta+h) - howfar(theta-h)) / (2h)
tstar = find_zero(howfarp, 45)
```

This graph shows the differences in the trajectories:

```figure
plot(x -> flight(x, 45), 0, howfar(45))  
plot!(x -> flight(x, tstar), 0, howfar(tstar))
```







## Potential issues

The higher-order methods are basically various derivative-free
versions of Newton's method (which has update step $x - f(x)/f'(x)$). For
example, Steffensen's method (`Order2()`) essentially replaces $f'(x)$ with
$(f(x + f(x)) - f(x))/f(x)$. This is a forward-difference
approximation to the derivative with "$h$" being $f(x)$, which
presumably is close to $0$ already. The methods with higher order
combine this with different secant line approaches that minimize the
number of function calls. These higher-order methods can be
susceptible to some of the usual issues found with Newton's method:
poor initial guess, small first derivative, or large second derivative
near the zero.

When the first derivative is near $0$, the value of the next step can
be quite different, as the next step generally tracks the intersection point of
the tangent line. We see that starting at a $\pi/2$ causes this search
to be problematic:

```
find_zero(sin, pi/2, Order1())
```

(Whereas, starting at `pi/2 + 0.3`--where the slope of the tangent
is sufficiently close to point towards $\pi$--will find convergence at $\pi$.)

For a classic example where a large second derivative is
the issue, we have $f(x) = x^{1/3}$:

```
f(x) = cbrt(x)
x = find_zero(f, 1, Order2())	# all of 2, 5, 8, and 16 fail or diverge towards infinity
```

However, the default finds the root here, as a bracket is identified:

```
x = find_zero(f, 1)
x,  f(x)
```

Finally, for many functions, all of these methods need a good initial
guess. For example, the polynomial function $f(x) = x^5 - x - 1$ has
its one zero near $1.16$. If we start far from it, convergence may
happen, but it isn't guaranteed:

```
f(x) = x^5 - x - 1
x0 = 0.1
find_zero(f, x0) 
```

Whereas, 

```
find_zero(f, x0, Order2())
```

A graph shows the issue. We have overlayed 15 steps of Newton's
method, the other algorithms being somewhat similar:

```figure,nocode
using ForwardDiff
D(f) = x -> ForwardDiff.derivative(f,float(x))
xs = [x0]
n = 15
for i in 1:(n-1) push!(xs, xs[end] - f(xs[end])/D(f)(xs[end])) end
ys = [zeros(Float64,n)';map(f, xs)'][1:2n]
xs = xs[repeat(collect(1:n), inner=[2], outer=[1])]
plot(f, -1.25, 1.5, linewidth=3, legend=false)
plot!(zero, -1.25, 1.5, linewidth=3)
plot!(xs, ys)
```

Though 15 steps are shown, only a few are discernible, as the
function's relative maximum causes a trap for this algorithm. Starting
to the right of the relative minimum--nearer the zero--would avoid
this trap. The default method employs a trick to bounce out of such
traps, though it doesn't always work.

###  Tolerances

Mathematically solving for a zero of a nonlinear function may be
impossible, so numeric methods are utilized. However, using floating
point numbers to approximate the real numbers leads to some nuances.


For example, consider the polynomial $f(x) = (3x-1)^5$ with one zero
at $1/3$ and its equivalent expression $f1(x) = -1 + 15\cdot x - 90\cdot
x^2 + 270\cdot x^3 - 405\cdot x^4 + 243\cdot x^5$. Mathematically
these are the same, however not so when evaluated in floating
point. Here we look at the 21 floating point numbers near $1/3$:

```
f(x) = (3x-1)^5
f1(x) =  -1 + 15*x - 90*x^2 + 270*x^3 - 405*x^4 + 243*x^5
ns = [1/3];
u=1/3; for i in 1:10 (u=nextfloat(u);push!(ns, u)) end
u=1/3; for i in 1:10 (u=prevfloat(u);push!(ns, u)) end
sort!(ns)

maximum(abs.(f.(ns) - f1.(ns)))
```

We see the function values are close for each point, as the maximum difference
is like $10^{15}$. This is roughly as expected, where even one
addition may introduce a relative error as big as $2\cdot 10^{-16}$ and here
there are several such.

Generally this is not even a thought, as the differences are generally
negligible, but when we want to identify if a value is zero, these
small differences might matter. Here we look at the signs of the
function values:


```
fs = sign.(f.(ns))
f1s = sign.(f1.(ns))
[ns-1/3 fs f1s]
```

Parsing this shows a few surprises. First, there are two zeros of
`f(x)` identified--not just one, as expected mathematically--the
floating point value of `1/3` and the next largest floating point
number. For `f1(x)` there is only one zero, but it isn't the floating
point value for `1/3` but rather 10 floating point numbers
away. Further, there are 5 sign changes of the function values. There
is no guarantee that a zero will be present, but for a mathematical
function that changes sign, there will be at least one sign change.

With this in mind, an "exact" zero of `f` would be either where
`iszero(f(x))` is true *or* where the function has a sign change at an
adjacent floating point value (either `f(x)*f(prevfloat(x))<0` or
`f(x)*f(nextfloat(x)) < 0`).

As mentioned, the default `Bisection()` method of `find_zero`
identifies such zeros for `f` provided an initial bracketing interval
is specified when `Float64` numbers are used.  However, if a
mathematical function does not cross the $x$ axis at a zero, then
there is no guarantee the floating point values will satisfy either of
these conditions.


Now consider the function `f(x) = exp(x)-x^4`. The value`x=8.613169456441398` is a zero in this sense, as there is a change of sign:

```
f(x) = exp(x) - x^4
F(x) = sign(f(x))
x=8.613169456441398
F(prevfloat(x)), F(x), F(nextfloat(x))
```

However, the value of `f(x)` is not as small as one might initially
expect for a zero:

```
f(x)
```

The value `x` is an approximation to the actual mathematical zero,
call it $x$. There is a difference between $f(x)$ (the mathematical answer) and `f(x)` (the floating point answer). Roughly speaking we expect `f(x)` to be about $f(x) + f'(x)\cdot \delta$, where $\delta$ is the difference between `x` and $x$. This will be on the scale of `abs(x) * eps()`, so all told we expect an answer to be in the range of $0$ plus or minus this value:

```
fp(x) = exp(x) - 4x^3 # the derivative
fp(x) * abs(x) * eps()
```

which is about what we see.

Bisection can be a slower method than others. For floating point values, `Bisection()` takes no more than 64 steps, but other methods may be able to converge to a zero in 4-5 steps (assuming good starting values are specified).

When fewer function calls are desirable, then checking for an
*approximate* zero may be preferred over assessing if a sign change
occurs, as generally that will take two additional function calls per
step. Besides, a sign change isn't guaranteed for all zeros. An approximate zero would be one where $f(x) \approx 0$.

By the above, we see that we must consider an appropriate
tolerance. The first example shows differences in floating point
evaluations from the mathematical ones might introduce errors on the
scale of `eps` regardless of the size of `x`. As seen in the second
example, the difference between the floating point approximation to
the zero and the zero introduces a potential error *proportional* to
the size of `x`. So a tolerance might consider both types of
errors. An absolute tolerance is used as well as a relative tolerance,
so a check might look like:

```verbatim
abs(f(x)) < max(atol, abs(x) * rtol)
```

This is different from `Julia`'s `isapprox(f(x), 0.0)`, as that would use `abs(f(x))` as the multiplier, which renders a relative tolerance useless for this question.

One issue with relative tolerances is that for functions with
sublinear growth, extremely large values will be considered zeros.
Returning to an earlier example, we have a misidentified zero:

```
find_zero(cbrt, 1, Order8())
```

For `Order8`, the algorithm rapidly marches off towards infinity so the relative
tolerance $\approx |x| \cdot \epsilon$ used to check if $f(x) \approx
0$ is large compared to the far-from zero $f(x)$.

Either the users must be educated about this possibility, or the
relative tolerance should be set to $0$. In that case, the absolute
tolerance must be relatively generous.  A conservative choice of
absolute tolerance might be `sqrt(eps())`, or about `1e-8`,
essentially the one made in SciPy.

This is not the choice made in `Roots`. The fact that bisection can
produce zeros as exact as possible, and the fact that the error in
function evaluation, $f'(x)|x|\epsilon$, is not typically on the scale
of `1e-8`, leads to a desire for smaller residual values of `f(x)` when possible.

In `Roots`, the faster algorithms use a check on both the size of
`f(xn)` and the size of the difference between the last two `xn` values. The check on `f(xn)`
is done with a tight tolerance, as is the check on $x_n \approx
x_{n-1}$. If the function values get close to zero, an
approximate zero is declared. Further, if the $x$ values get close to each other
*and* the function value is close to zero with a *relaxed* tolerance,
then an approximate zero is declared. In practice this seems to work
reasonably well. The relaxed tolerance uses the cube root of the
absolute and relative tolerances.







## Searching for all zeros in an interval

The methods described above are used to identify one of possibly
several zeros.  The `find_zeros` function searches the interval $(a,b)$
for all zeros of a function $f$. It is straightforward to use:

```
f(x) = exp(x) - x^4
a, b = -10, 10
zs = find_zeros(f, a, b)
```

The  search interval, $(a,b)$, is specified through two arguments. It is
assumed that neither endpoint is a zero. Here we see the result of the
search graphically:

```
plot(f, a, b)
scatter!(zs, f.(zs))
```

We can identify points where the first and second derivative is
zero. We use `D` from above:


```
f(x) = cos(x) + cos(2x)
a, b = -10, 10
cps = find_zeros(D(f), a, b)
ips = find_zeros(D(f,2), a, b)
plot(f, a, b)
scatter!(cps, f.(cps))
scatter!(ips, f.(ips), markercolor = :yellow)
```

The `find_zeros` algorithm will use bisection when a bracket is
identified. This method will identify jumps, so areas where the
derivative changes sign (and not necessarily a zero of the derivative)
will typically be identified:

```
f(x) = abs(cbrt(x^2-1))
a, b = -5, 5
cps = find_zeros(D(f), a, b)
plot(f, a, b)
scatter!(cps, f.(cps))
```

In this example, the derivative has vertical asymptotes at $x=1$ and
$x=-1$ so is not continuous there. The bisection method identifies the
zero crossing, not a zero.




The search for all zeros in an interval is confounded by a few things:

* too many zeros in the interval $(a,b)$
* nearby zeros

The algorithm is adaptive, so that it can succeed when there are many
zeros, but it may be necessary to increase `no_pts` from the default
of 12, at the cost of possibly taking longer for the search.

Here the algorithm identifies all the zeros, despite there being several:

```
f(x) = cos(x)^2 + cos(x^2)
a, b = 0, 10
rts = find_zeros(f, a, b)
plot(f, a, b)
scatter!(rts, f.(rts))
```

For nearby zeros, the algorithm does pretty well, though it isn't
perfect.

Here we see for $f(x) = \sin(1/x)$--with infinitely many zeros around
$0$--it finds many:

```
f(x) = iszero(x) ? NaN : sin(1/x)  # avoid sin(Inf) error
rts = find_zeros(f, -1, 1)  # 88 zeros identified
plot(f, -1, 1)
scatter!(rts, f.(rts))
```

The function, $f(x) = (x-0.5)^3 \cdot (x-0.499)^3$, looks *too* much like
$g(x) = x^6$ to `find_zeros` for success, as the two zeros are very nearby:

```
f(x) =  (x-0.5)^3 * (x-0.499)^3
find_zeros(f, 0, 1)
```

The issue here isn't *just* that the algorithm can't identify zeros
within $0.001$ of each other, but that the high power makes many
nearby values approximately zero.

The algorithm will have success when the powers are smaller

```
f(x) =  (x-0.5)^2 * (x-0.499)^2
find_zeros(f, 0, 1)
```

It can have success for closer pairs of zeros:

```
f(x) = (x-0.5) * (x - 0.49999)
find_zeros(f, 0, 1)
```

Combinations of large (even) multiplicity zeros or very nearby
zeros, can lead to misidentification.
