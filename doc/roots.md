# The Roots package

The `Roots` package contains simple routines for finding zeros of
continuous scalar functions of a single real variable.  A zero of $f$
is a value $c$ where $f(c) = 0$.  The basic interface is through the
function `find_zero`, which through multiple dispatch can handle many different cases.

We will use `Plots` for plotting.

```
using Roots  
using Plots
pyplot()
```

## Bracketing

For a function $f: R \rightarrow R$ a *bracket* is a pair $ a < b $ for
which $f(a) \cdot f(b) < 0$. That is they have different signs. If $f$
is a continuous function this ensures there to be a zero (a $c$ with
$f(c) = 0$) in the interval $[a,b]$,
otherwise, if $f$ is only piecewise continuous, there must be a point
$c$ in $[a,b]$ with the left limit and right limit at $c$ having
different signs (or $0$). Such values can be found, up to floating point
roundoff.

That is, given `f(a) * f(b) < 0`, a value `c` with `a < c < b` can be
found where either `f(c) == 0.0` or at least `f(prevfloat(c)) *
f(nextfloat(c)) <= 0`.

To illustrate, consider the function $f(x) = \cos(x) - x$. From the
graph we see readily that $[0,1]$ is a bracket (which we emphasize
with an overlay):

```figure
f(x) = cos(x) - x
plot(f, -2, 2)
plot!([0,1], [0,0], linewidth=2)
```

We use a vector or tuple to specify the initial condition for `Bisection`:

```
x = find_zero(f, (0, 1), Bisection())    # alternatively fzero(f, [0, 1])
x, f(x)
```

For this function we see that `f(x) == 0.0`.

----

Next consider $f(x) = \sin(x)$. A known root is $\pi$. Trignometry
tells us that $[\pi/2, 3\pi/2]$ will be a bracket. Here `Bisection()`
is not specified, as it will be the default when the initial value is
specified as pair of numbers:

```
f(x) = sin(x)
x = find_zero(f, (pi/2, 3pi/2))
x, f(x)
```

This value of `x` does not produce `f(x) == 0.0`, however, it is as close as can be:

```
f(prevfloat(x)) * f(x) < 0.0 || f(x) * f(nextfloat(x)) < 0.0
```

That is, at `x` the function is changing sign.

From a mathematical perspective, a zero is guaranteed for a
*continuous* function. However, the computer algorithm doesn't assume
continuity, it just looks for changes of sign. As such, the algorithm
will  identify discontinuities, not just zeros. For example:

```
find_zero(x -> 1/x, (-1, 1))
```


The endpoints can even be infinite:

```
find_zero(x -> Inf*sign(x), (-Inf, Inf))  # Float64 only
```


The basic algorithm used for bracketing when the values are simple
floating point values is the bisection method. For big float values, an
algorithm due to Alefeld, Potra, and Shi is used.

```
find_zero(sin, (big(3), big(4)))    # uses a different algorithm then for (3,4)
```

## Using an initial guess

Bracketing methods have guaranteed convergence, but in general require
many more function calls than are needed to produce an answer.  If a
good initial guess is known, then the `find_zero` function provides an
interface to some different iterative algorithms that are more
efficient. Unlike bracketing methods, these algorithms may not
converge to the desired root if the initial guess is not well chosen.

The default algorithm is modeled after an algorithm used for
[HP-34 calculators](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf). This
algorithm is designed to be more forgiving of the quality of the
initial guess at the cost of possibly performing many more steps. In
many cases it satisfies the criteria for a bracketing solution, as it
will use bracketing if within the algorithm a bracket is identified.

For example, the answer to our initial problem is near 1. Given this,
we can find the zero with:

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

The default call to `fzero` uses a first order method and
then possibly bracketing, which involves potentially many more function
calls. Though specifying a initial value is more convenient than a
bracket, there may be times where a more efficient algorithm is sought.
For such, a higher-order method might be better
suited. There are algorithms `Order1` (secant method), `Order2`
([Steffensen](http://en.wikipedia.org/wiki/Steffensen's_method)), `Order5`,
`Order8`, and `Order16`. The order 2 method is generally more efficient, but is more
sensitive to the initial guess than, say, the order 8 method. These
algorithms are accessed by specifying the method after the initial point:

```
f(x) = 2x - exp(-x)
x = find_zero(f, 1, Order2())      # also fzero(f, 1, order=2)
x, f(x)
```

```
f(x) = (x + 3) * (x - 1)^2
x = find_zero(f, -2, Order5())
x, f(x)
```

```
x = find_zero(f, 2, Order8())
x, f(x)
```

The latter shows that zeros need not be simple zeros (i.e. $f'(x) = 0$,
if defined) to be found. (Though non-simple zeros may take many more
steps to converge.)

To investigate the algorithm and its convergence, the argument
`verbose=true` may be specified.


For some functions, adjusting the default tolerances may be necessary
to achieve convergence. These include `atol` and `rtol`, which are
used to check if $f(x_n) \approx 0$;
`xatol`, `xrtol`, to check if $x_n \approx x_{n-1}$; and `maxevals` and `maxfnevals` to limit the
number of steps in the algorithm or function calls.


The higher-order methods are basically various derivative-free
versions of Newton's method (which has update step $x - f(x)/f'(x)$). For
example, Steffensen's method is essentially replacing $f'(x)$ with
$(f(x + f(x)) - f(x))/f(x)$. This is a forward-difference
approximation to the derivative with "$h$" being $f(x)$, which
presumably is close to $0$ already. The methods with higher order
combine this with different secant line approaches that minimize the
number of function calls. These higher-order methods can be
susceptible to some of the usual issues found with Newton's method:
poor initial guess, small first derivative, or large second derivative
near the zero.

----

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

Order 8 illustrates that sometimes the stopping rules can be misleading and
checking the returned value is always a good idea:

```
find_zero(f, 1, Order8())
```

The algorithm rapidly marches off towards infinity so the relative
tolerance $\approx |x| \cdot
\epsilon$ is large compared to the far-from zero $f(x)$.

----

This example illustrates that the default `find_zero`
call is more forgiving to an initial guess. The devilish function
defined below comes from a [test
suite](http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html)
of difficult functions. The default method finds the zero starting at 0:

```figure
f(x) = cos(100*x)-4*erf(30*x-10)
plot(f, -2, 2)
```

```
find_zero(f, 0)
```

Whereas, with higher order methods fail. For example,

```
find_zero(f, 0, Order8())
```

Basically the high order oscillation can send the proxy tangent line
off in nearly random directions. The default method can be fooled here
too.

----

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
D(f) = x -> ForwardDiff.derivative(f,x)
xs = [x0]
n = 15
for i in 1:(n-1) push!(xs, xs[end] - f(xs[end])/D(f)(xs[end])) end
ys = [zeros(Float64,n)';map(f, xs)'][1:2n]
xs = xs[repeat(collect(1:n), inner=[2], outer=[1])]
plot(f, -1.25, 1.5, linewidth=3, legend=false)
plot!(zero, -1.25, 1.5, linewidth=3)
plot!(xs, ys)
```

Though 15 steps are shown, only a few are discernible, as the function's relative maximum
causes a trap for this algorithm. Starting to the right of the
relative minimum -- nearer the zero -- would avoid this trap. The default
method employs a trick to bounce out of such traps, though it doesn't
always work.

## Finding more than one zero

The bracketing methods suggest a simple algorithm to recover multiple
zeros: partition an interval into many small sub-intervals. For those
that bracket a root find the root. This is essentially implemented with `find_zeros(f,
a, b)`. The algorithm has problems with non-simple zeros (in particular ones
that don't change sign at the zero) and zeros which bunch together. Simple usage is often succesful
enough, but a graph should be used to assess if all the zeros are found:

```
find_zeros(x -> exp(x) - x^4, -10, 10)
```



## Classical methods

The package provides some classical methods for root finding:
`newton`, `halley`, and `secant_method`. We can see how each works on
a problem studied by Newton himself. Newton's method uses the function
and its derivative:

```
f(x) = x^3 - 2x - 5
fp(x) = 3x^2 - 2
x = newton(f, fp, 2)
x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
```

To see the algorithm in progress, the argument `verbose=true` may be specified. 

The secant method needs two starting points, here we start with 2 and 3:

```
x = secant_method(f, 2,3)
x, f(x), f(prevfloat(x)) * f(nextfloat(x))
```

Halley's method has cubic convergence, as compared to Newton's
quadratic convergence. It uses the second derivative as well:

```
fpp(x) = 6x
x = halley(f, fp, fpp, 2)
x, f(x), f(prevfloat(x)) * f(nextfloat(x))
```


For many functions, the derivatives can be computed automatically. The
`ForwardDiff` package provides a means. Here we define an operator `D`
to compute a derivative:

```
using FowardDiff
D(f) = x -> ForwardDiff.derivative(f, x)
D(f, n) = n > 1 ? D(D(f),n-1) : D(f)
```


```
newton(f, D(f), 2)    
```

Or for Halley's method

```
halley(f, D(f), D(f,2), 2)  
```


## Finding critical points

The `D` function makes it straightforward to find critical points
(where the derivative is $0$ or undefined). For example, the critical
point of the function $f(x) = 1/x^2 + x^3, x > 0$ near $1.0$ can be found with:

```
f(x) = 1/x^2 + x^3
find_zero(D(f), 1)
```

For more complicated expressions, `D` will not work. In this example,
we have a function $f(x, \theta)$ that models the flight of an arrow
on a windy day:

```
function flight(x, theta)
 	 k = 1/2
	 a = 200*cosd(theta)
	 b = 32/k
	 tand(theta)*x + (b/a)*x - b*log(a/(a-x))
end
```



The total distance flown is when `flight(x) == 0.0` for some `x > 0`:
This can be solved for different `theta` with `fzero`. In the
following, we note that `log(a/(a-x))` will have an asymptote at `a`,
so we start our search at `a-5`:

```
function howfar(theta)
	 a = 200*cosd(theta)
	 find_zero(x -> flight(x, theta), a-5)
end
```	 

To see the trajectory if shot at 45 degrees, we have:


```figure
theta = 45
plot(x -> flight(x,  theta), 0, howfar(theta))
```


To maximize the range we solve for the lone critical point of `howfar`
within the range. The derivative can not be taken automatically with
`D`. So,  here we use a central-difference approximation and start the
search at 45 degrees, the angle which maximizes the trajectory on a
non-windy day:

```
h = 1e-5
howfarp(theta) = (howfar(theta+h) - howfar(theta-h)) / (2h)
tstar = find_zero(howfarp, 45)
```

This graph shows the differences in the trajectories:

```figure
plot(x -> flight(x, tstar), 0, howfar(tstar))
plot!(x -> flight(x, 45), 0, howfar(45))
```



# Use with other number types


The `Unitful` package provides a means to attach units to numeric
values.

For example, a projectile motion with $v_0=10$ and $x_0=16$ could be
represented with:

```
using Unitful
s = u"s"; m = u"m"
g = 9.8*m/s^2
v0 = 10m/s
y0 = 16m
y(t) = -g*t^2 + v0*t + y0
```

This motion starts at a height of 16 meters and has an initial
velocity of 10 meters per second.

The time of touching the ground is found with:

```
a = find_zero(y, 1s, Order2())
a
```

Automatic derivatives don't propogate through `Unitful`, so we define
the approximate derivative--paying attention to units--with:

```
Df(f, h=1e-6) = x -> (f(x + h*oneunit(x)) - f(x)) / (h*oneunit(x))
```

And then the fact the peak is the only local maximum, it can be found from:

```
find_zero(Df(y), a/2, Order2())
```

----

The `SymEngine` package provides symbolic values to `Julia`. Rather
than passing a function to `find_zero`, we can pass a symbolic expression:

```
using SymEngine
g, v0, y0 = 9.8, 10, 16
@vars t
yt = -g * t^2 + v0 * t + y0
```

```
a = find_zero(yt, 1, Order2())
a
```

And the peak is determined to be at:

```
find_zero(diff(yt, t), (0, a), Bisection())
```

(This also illustrates that symbolic values can be passed to describe
the `x`-axis values.)

Similarly, the `SymPy` package could be used in an identical way.
