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
is a continuous function this ensures there will be a zero
in the interval $[a,b]$.
Otherwise, if $f$ is only piecewise continuous, there must be a point
$c$ in $[a,b]$ with the left limit and right limit at $c$ having
different signs (or $0$).

Such values can be found, up to floating point
roundoff. That is, given `f(a) * f(b) < 0`, a value `c` with `a < c < b` can be
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
will  identify discontinuities, not just zeros. For example, here we
see the vertical asymptote identified:

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

By default, bisection will converge to machine tolerance. This may
provide more accurancy than desired. A tolerance may be specified to
terminate early, thereby utilizing fewer resources. For example, the following
uses 19 steps to reach accuracy to $10^{-6}$ (without specifying `xatol` it uses
52 steps):

```
rt =find_zero(sin, (3.0, 4), xatol=1e-6)
rt - pi
```


## Non-bracketing algorithms

Bracketing methods have guaranteed convergence, but in general
identifying a bracket can take  more effort and the algorithms may require
many more function calls than are needed to produce an answer.  If a
good initial guess is known, then the `find_zero` function provides an
interface to some different iterative algorithms that are usually more
efficient. Unlike bracketing methods, these algorithms may not
converge to the desired root if the initial guess is not well chosen.

The default algorithm is modeled after an algorithm used for
[HP-34 calculators](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf). This
algorithm is designed to be more forgiving of the quality of the
initial guess at the cost of possibly performing many more steps. If
the algorithm encounters a bracket, bisection will be used, so it may
use many function calls.

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

The default call to `fzero` uses a first order method and then
possibly bracketing, which involves potentially many more function
calls. Though specifying a initial value is more convenient than a
bracket, there may be times where a more efficient algorithm is
sought.  For such, a higher-order method might be better suited. There
are algorithms `Order1` (secant method), `Order2`
([Steffensen](http://en.wikipedia.org/wiki/Steffensen's_method)),
`Order5`, `Order8`, and `Order16`.  The order 1 or 2 methods are
generally quite efficient. The higher order ones may be more efficient
when higher precision is sought.  These algorithms are accessed by
specifying the method after the initial starting point:

```
f(x) = 2x - exp(-x)
x = find_zero(f, 1, Order1())      # also fzero(f, 1, order=1)
x, f(x)
```

The above makes $8$ function calls, to the $57$ made with `Order0`.

```
f(x) = (x + 3) * (x - 1)^2
x = find_zero(f, -2, Order2())
x, f(x)
```


```
x = find_zero(f, 2, Order8())
x, f(x)
```

The latter shows that zeros need not be simple zeros (i.e. $f'(x) =
0$, if defined) to be found. Generally speaking, non-simple zeros are
expected to take many more function calls, as the methods are no
longer super-linear. This is the case here where `Order2` uses $55$
function calls, `Order8` uses $41$, and `Order0` takes a comparable $42$.)

To investigate an algorithm and its convergence, the argument
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

Relative tolerances are important when answers could possible be quite
large. Just evaluating a floating point approximation to a zero can
produce a value on the size of $f'(x) x \epsilon$. If it is known that
large $x$ values are not sought, then setting the relative tolerance
to zero will avoid the issue of marching off towards infinity and
returning what appears to be a valid zero:

```
find_zero(f, 1, Order8(), rtol=0.0)
```

----

This example illustrates that the default `find_zero`
call is more forgiving to an initial guess. The devilish function
defined below comes from a [test
suite](http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html)
of difficult functions. The default method finds the zero starting at 0:

```figure
using SpecialFunctions
f(x) = cos(100*x)-4*erf(30*x-10)
plot(f, -2, 2)
```

```
find_zero(f, 0)
```

Whereas, with higher order methods fail. For example,

```
find_zero(f, 0, Order1())
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

Though 15 steps are shown, only a few are discernible, as the function's relative maximum
causes a trap for this algorithm. Starting to the right of the
relative minimum -- nearer the zero -- would avoid this trap. The default
method employs a trick to bounce out of such traps, though it doesn't
always work.


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
x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
```

To see the algorithm in progress, the argument `verbose=true` may be
specified.

Alternatively, `Roots.Newton()` can be specified as the method for `find_zero`. The
functions are specified as a tuple:

```
find_zero((f,fp), 2, Roots.Newton())
```

The secant method typically needs two starting points, here we start with 2 and 3:

```
x = Roots.secant_method(f, (2,3))
x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
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


For many functions, the derivatives can be computed automatically. The
`ForwardDiff` package provides a means. Here we define an operator `D`
to compute a derivative:

```
using FowardDiff
D(f) = x -> ForwardDiff.derivative(f, float(x))
D(f, n) = n > 1 ? D(D(f),n-1) : D(f)
```


```
Roots.newton(f, D(f), 2)    
```

Or for Halley's method

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
of finding a derivative can be empoloyed. In

this example, we have a function $f(x, \theta)$ that models the flight
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
plot(x -> flight(x, 45), 0, howfar(45))  
plot!(x -> flight(x, tstar), 0, howfar(tstar))
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

```
plot(yt, 0, a)
```



And the peak is determined to be at:

```
find_zero(diff(yt, t), (0, a), Bisection())
```

(This also illustrates that symbolic values can be passed to describe
the `x`-axis values.)

Similarly, the `SymPy` package could be used in an identical way.



## Searching for all zeros in an interval

The methods described above are used to identify one of possible
several zeros.  The `find_zeros` function searches interval $(a,b)$
for all zeros of a function $f$. It is straightforward to use:

```
f(x) = exp(x) - x^4
a, b = -10, 10
zs = find_zeros(f, a, b)
```

The  search interval $(a,b)$ is specified through two arguments. It is
assumed that neither endpoint is a zero. Here we see the result of the
search graphically:

```
plot(f, a, b)
scatter!(zs, f.(zs))
```

We can identify points where the first and second derivative is
zero. We use `D` from above


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

----

Leveraging the zero-finding tools available in `Roots`, the task of finding all zeros in an interval $(a,b)$ can be approached by a naive method:

* split the interval into $n$ points $a=x_1 < x_2 < \dots < x_n=b$.

* There are $n-1$ subintervals made from adjacent pairs. Each subinterval, $(x_i, x_{i+1})$, is considered: if it is a bracketing interval, *bisection* is used to find a guaranteed zero; if not, `find_zero` is used to search for a zero in $(x_i, x_{i+1})$.

* collect the zeros found

Bisection is guaranteed to find a zero, so the cost of many function calls is worth paying. But calls to  `find_zero` for non-bracketing intervals have no guarantee and may be foolish. Though the cost (in function calls) is not generally as high as bisection, it is best to avoid excessive calls. As such, the `find_zeros` implementation oversamples the subintervals checking for bracketing intervals. This reduces the amount of resources `find_zeros` needs.
Still, when a function is "expensive" to compute the approach can use too many function calls. (The above two examples use between 1000 and 6000 function calls.)


This approach can be useful, but certain scenarios will be problematic:


* Bisection will only find one answer in a bracketing interval, but a continuous function can have one, three, five, $\dots$ zeros.

* A non bracketing interval may contain  zeros (even an even number crossing the $x$ axis) but `find_zero` may fail to find any, as convergence can be very dependent on the initial starting point.

* Nearby zeros can cause issues. For example with $f(x) = (x-0.5)^3 \cdot (x - 0.499)^3$ the function looks a lot like $f(x) = (x-0.5)^6$ and the algorithm might just find one zero. (In fact, `find_zeros` fails on this problem, but will have luck with $f(x) = (x-0.5)^2 \cdot (x - 0.499)^2$.)

* A highly oscillatory function may have many more zeros than the initally specified number of subintervals


For example, consider the cosine curve with period $1$:

```
f(x) = cos(2*pi*x)  # is 1.0 at 0,1,2,3,....
```


This function has a zero at $1/4, 3/4, 5/4, \dots$; altogether there
are $32$ zeros in $(0,16)$. However, if we break the interval $(0,16)$
into 4 subintervals, and oversample each using $4$ subintervals, we
will be looking at the function at $0,1,2,\dots, 16$ where it is
always $1$. So despite the many oscillations, no bracketing intervals
will be found. Furthermore if we start the derivative free search at
the midpoint of the `no_pts-1` intervals, it will start at $2$, $6$,
$10$, and $14$. For these values, again the function is always at its
maximum of $1$. Relative extrema are especially poor choices of
starting point, as the tangent lines have $0$ slope. So we expect the
algorithm might have issues, and it does:

```
find_zeros(f, 0, 16, no_pts=5, k=4)  # no pts - 1 = no subintervals
```

We see the algorithm missed many of the answers. The default values won't miss, but a similar, systematically bad example could be constructed for them. So users must be wary and may need to modify the interval selection as above (where the number of initial subintervals is specified through `no_pts`, it being `no_pts-1`, and the number of oversampled intervals considered for bracketing is specified by $k$.

Consider this choice, with just one additional subinterval for the
bisection check. Now all $32$ answers are found:

```
pts = find_zeros(f, 0, 16, no_pts=5, k=5)
plot(f, 0, 16)
scatter!(pts, f.(pts))
```



An analysis shows that now the initial points considered for brackets
are $0$ through $16$ with a step size of $0.8$ (and not $1.0$). This
step size avoids the aliasing and the algorithm is successful.

Generally though, if it is suspected that the number of zeros
identified is too few, increasing `no_pts` from its default of 11 will
increase the chance of finding more zeros. 

The above uses only 20 intervals initially, yet $32$ answers are
found. How? When zeros are found, after stepping away from the zeros,
the algorithm will recursively try again on the subintervals between
the zeros found, though using a reduced number of points. This allows
reasonably nearby points to be identified in a subsequent pass and can
find additional zeros when at least one is found using the initial
intervals.

Of course, this recursive calling doesn't go for ever. For example,
the function $f(x) = \sin(1/x)$ has inifinitely many zeros near
$0$. The algorithm identifies many, but could never find all:

```
f(x) = iszero(x) ? NaN : sin(1/x)  # redefine to avoid an error with `sin(Inf)`
rts = find_zeros(f, -1, 1)
length(rts)
```

What is found are zeros (with maximum function value around `1e-13`):

```
maximum(abs.(f.(rts)))
```

The resolution does get fairly small:

```
minimum(diff(rts))  # smallest gap 
```


The closest nearby roots are on the scale of $10^{-4}$, which is greater than the resolution of the graph, when zoomed in:

```
plot(f, -1, 1)
scatter!(rts, f.(rts))
```
