# The Roots package

The `Roots` package contains simple routines for finding roots of continuous scalar functions of a single real variable.  The basic interface is through the function `fzero`, which through multiple dispatch can handle many different cases.

We will use these pacakges

```
using Gadfly			# for plotting
using Roots
using Polynomial
```

## Bracketing

For a function $f: R \rightarrow R$ a bracket is a pair $a<b$ for which $f(a)\cdot f(b) < 0$. That is they have different signs. If $f$ is continuous this forces there to be a zero on the interval $[a,b]$, otherwise, if $f$ is only piecewise continuous, there must be a point $c$ in $[a,b]$ with the left limit and right limit at $c$ having different signs. These values can be found, up to floating point roundoff.

That is, a value $a < c < b$ can be found with either `f(c) == 0.0` or `f(prevfloat(c)) * f(nextfloat(c)) <= 0`.

To illustrate, consider the function $f(x) = \cos(x) - x$. From the graph we see readily that $[0,1]$ is a bracket:

```
f(x) = cos(x) - x
plot(f, -2,2)
```

The basic function call specifies a bracket using vector notation:

```
x = fzero(f, [0, 1])
[x, f(x)]
```

For that function `f(x) == 0.0`. Next consider $f(x) = \sin(x)$. A known root is $\pi$. Basic trignometry tells us that $[\pi/2, 3\pi2]$ will be a bracket:

```
f(x) = sin(x)
x = fzero(f, [pi/2, 3pi/2])
[x, f(x)]
```

This value of `x` does not produce `f(x) == 0.0`, however, it is as close as can be:

```
f(prevfloat(x)) * f(nextfloat(x)) < 0.0
```

That is at `x` the function is changing sign.

(The latter guarantee should actually just be either `f(x)*f(prevfloat(x)) < 0` or `f(x) * f(nextfloat(x)) < 0`, but that is a bit much to write, and only a very deviously defined function could be an exception.)

The basic algorithm used for bracketing when the values are simple floating point values is the bisection method. Though there are algorithms that mathematically should converge faster (and one is used for the case where `BigFloat` numbers are used) by exploiting floating point computations this algorithm uses fewer function calls and runs faster.

## Using an initial guess

If a bracket is not known, but a good initial guess is, the `fzero` provides an interface to some different algorithms. The basic algorithm is modeled after an algorithm used for [HP-34 calculators](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf). This algorithm is much more robust to the quality of the initial guess and does not rely on tolerances for a stopping rule. In many cases it satisfies the criteria for bracketing.

For example, we have:

```
f(x) = cos(x) - x
x = fzero(f , 1)
[x, f(x)]
```

And 

```
f(x) = x^3 - 2x - 5
x = fzero(f, 2)
[x, f(x), f(prevfloat(x)) * f(nextfloat(x))]
```

For even more precision, `BigFloat` numbers can be used

```
x = fzero(sin, big(3))
[x, f(x), x - pi]
```

### Higher order methods

The default call to `fzero` uses a second order method at best and then bracketing, which involves potentially many more function calls. For some functions, a higher-order method might be better suited. There are algorithms of order 1 (secant method), 2 ([Steffensen](http://en.wikipedia.org/wiki/Steffensen's_method)), 5, 8, and 16. The order 2 method is generally more efficient, but is more sensitive to the initial guess than, say, the order 8 method. These algorithms are accessed by specifying a value for the `order` argument:

```
f(x) = 2x - exp (-x)
x = fzero(f, 1, order=2)
[x, f(x)]
```

```
f(x) = (x + 3) * (x - 1)^2
x = fzero(f, -2, order=5)
(x, f(x))
```

```
x = fzero(f, 2, order=8)
(x, f(x))
```

The latter shows that zeros need not be simple zeros (i.e. $f'(x) = 0$ if defined) to be found. For the higher-order methods, there is a tolerance that can be specified so that a value is returned as a zero if `abs(f(x)) < tol`. The default method for `fzero` uses a very strict tolerance for this, otherwise defaulting to an error that at times might be very close to the actual zero. For this problem it finds the exact value:

```
x = fzero(f, 2)
(x, f(x))
```

But not for a similar problem:

```
fzero(x -> x^4, 1)
```

(Though the answer is basically on track, the algorithm takes too long to improve itself to the very stringent range set. For problems where a bracket is found, this dithering won't happen.)


The higher-order methods are basically various derivative-free versions of Newtons method which has update step $x - f(x)/f'(x)$. For example, Steffensen's method is essentially replacing $f'(x)$ with $(f(x + f(x)) - f(x))/f(x)$. This is just a forward-difference approximation to the derivative with "$h$" being $f(x)$, which presumably is close to $0$ already. The methods with higher order combine this with different secant line approaches that minimize the number of function calls. The default method uses a combination of Steffensen's method with modifications, a quadratic fit, and, if possible, a bracketing approach. It may need many more function calls than the higher-order methods. These higher-order methods can be susceptible to some of the usual issues found with Newton's method: poor initial guess, small first derivative, or large second derivative at the zero.


For a classic example where basically the large second derivative is the issue, we have $f(x) = x^{1/3}$:

```
f(x) = cbrt(x)
x = fzero(f, 1, order=8)	# all of 2, 5, 8, and 16 fail
(x, f(x))
```

However, the default finds the root here

```
x = fzero(f, 1)
(x, f(x))
```

Finally, we show another example illustrating that the default `fzero` call is more forgiving to an initial guess. The devilish function defined below comes from a [test suite](http://people.sc.fsu.edu/~jburkardt/cpp_src/test_zero/test_zero.html) of difficult functions. The default method finds the zero:

```
f(x) = cos(100*x)-4*erf(30*x-10)
plot(f, -2, 2)
```

```
fzero(f, 1)
```

Whereas, with `order=8` an error is raised.

```
fzero(f, 1, order=8)
```

In this example, some other orders work. Basically the high order oscillation can send the proxy tangent line off in nearly random directions.

## Polynomials

The `Polynomial` package provides a type for working with polynomial functions that allows many typical polynomial operations to be defined. In this context, the `roots` function is used to find the roots of a polynomial.

The `fzero` function will also find the roots. The algorithm does a better job when there are multiple roots, as it implements an algorithm that first identifies the multiplicity structure of the roots, and then tries to improve these values.

For example, 

```
x = Poly([1.0,0])			# "1x + 0"
roots((x-1)*(x-2)*(x-3))
```

The same is done with `fzero`:

```
fzero((x-1)*(x-2)*(x-3))
```

The `roots` function degrades as there are multiplicities:

```
p = (x-1)^2*(x-2)^3*(x-3)^4
roots(p)
```

Whereas,

```
fzero(p)
```

The difference gets dramatic when the degree is quite large (which is more of a mathematical possibility perhaps than a practical concern):

```
p = (x-1)^20 * (x-2)^5
r = roots(p)
plot(x = real(r), y=imag(r))
```

But,

```
fzero(p)
```


The `roots` function does not produce the accuracy that can be found using one of the algorithms of `fzero` for functions. For example

```
p = (x - sqrt(5))^3
r = roots(p)
r[1] - sqrt(5)
```

The `fzero` method can improve this, but if only one value is sought, a starting point can be given and the accuracy improved:

```
x = fzero(p, 2)
(x, polyval(p, x))
```



## Classical methods

The package provides some classical methods for root finding: `newton`, `halley`, and `secant_method`. We can see how each works on a problem studied by Newton himself. Newton's method uses  the function and its derivative:

```
f(x) = x^3 - 2x - 5
fp(x) = 3x^2 - 2
x = newton(f, fp, 2)
[x, f(x), f(prevfloat(x)) * f(nextfloat(x))]
```

To see the algorithm in progress, the argument `verbose=true` may be specified. In this case only 4 steps are needed:

```
newton(f, fp, 2, verbose=true)
```

The secant method needs two starting points, here we start with 2 and 3:

```
x = secant_method(f, 2,3, verbose=true)
[x, f(x), f(prevfloat(x)) * f(nextfloat(x))]
```

Halley's method has cubic convergence, as compared to Newton's quadratic convergence. It uses the second derivative as well:

```
fpp(x) = 6x
x = halley(f, fp, fpp, 2, verbose=true)
[x, f(x), f(prevfloat(x)) * f(nextfloat(x))]
```


For many function, the derivatives can be computed automatically. The `PowerSeries` function provides a means. This package wraps the process into an operator, `D` which returns the derivative of a function `f`:

```
newton(f, D(f), 2)
```

Or for Halley's method

```
halley(f, D(f), D(f,2), 2)
```

(The operator `D2(f)` is a convenience for `D(f,2)`.) Specifying the derivative(s) can be skipped, the functions will default to the above calls.

## Finding critical points

The `D` function makes it straightforward to find critical points (where the derivative is $0$ or undefined). For example, the critical point of the function $f(x) = 1/x^2 + x^3, x > 0$ can be found with:

```
f(x) = 1/x^2 + x^3
fzero(D(f), 1)
```

For more complicated expressions, `D` will not work. In this example, we have a function $f(x, theta)$ that models the flight of an arrow on a windy day:

```
function flight(x, theta)
 	 k = 1/2
	 a = 200*cosd(theta)
	 b = 32/k
	 tand(theta)*x + (b/a)*x - b*log(a/(a-x))
end
```

The total distance flown is when `flight(x) == 0.0` for some `x > 0`: This can be solved for different `theta` with `fzero`. In the following, we note that `log(a/(a-x))` will have an asymptote at `a`, so we start our search at `a-1`:

```
function howfar(theta)
	 a = 200*cosd(theta)
	 fzero(x -> flight(x, theta), a-1)
end
```	 

To maximize this function, the derivative can not be taken by `D`. Here we use a central-difference approximation and start the search at 45 degrees, the angle which maximizes the trajectory on a non-windy day:

```
h = 1e-5
howfarp(theta) = (howfar(theta+h) - howfar(theta-h)) / (2h)
fzero(howfarp, 45)
```
