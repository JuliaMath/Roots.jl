# An overview of `Roots`

The `Roots` package contains simple routines for finding zeros of
continuous scalar functions of a single real variable.  A zero of $f$
is a value $c$ where $f(c) = 0$.  The basic interface is through the
function `find_zero`, which through multiple dispatch can handle many different cases.

The [NonlinearSolve](https://github.com/JuliaComputing/NonlinearSolve.jl) package provides an alternative.

In the following, we will use `Plots` for plotting and `ForwardDiff`
to take derivatives.

```jldoctest roots
julia> using Roots

julia> using Plots, ForwardDiff

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

Such values can be found, up to floating point
roundoff. That is, given `f(a) * f(b) < 0`, a value `c` with `a < c < b` can be
found where either `f(c) == 0.0` or  `f(prevfloat(c)) * f(c) < 0` or
`f(c) * f(nextfloat(c)) < 0`.

To illustrate, consider the function $f(x) = \cos(x) - x$. From a
graph we can see readily that $[0,1]$ is a bracket

The `Roots` package includes the bisection algorithm through
`find_zero`. We use a structure for which `extrema` returns `(a,b)`
with `a < b`, such as a vector or tuple, to specify the initial
condition and `Bisection()` to specify the algorithm:

```jldoctest roots
julia> f(x) = cos(x) - x
f (generic function with 1 method)

julia> x = find_zero(f, (0, 1), Bisection())    # alternatively fzero(f, [0, 1])
0.7390851332151607

julia> x, f(x)
(0.7390851332151607, 0.0)
```

For this function we see that `f(x)` is `0.0`.

----

Next consider $f(x) = \sin(x)$. A known zero is $\pi$. Trigonometry
tells us that $[\pi/2, 3\pi/2]$ will be a bracket. In this call `Bisection()`
is not specified, as it will be the default when the initial value is
specified as a pair of numbers:

```jldoctest roots
julia> f(x) = sin(x)
f (generic function with 1 method)

julia> x = find_zero(f, (pi/2, 3pi/2))
3.141592653589793

julia> x, f(x)
(3.141592653589793, 1.2246467991473532e-16)

```

This value of `x` does not exactly produce a zero, however, it is as close as can be:

```jldoctest roots
julia> f(prevfloat(x)) * f(x) < 0.0 || f(x) * f(nextfloat(x)) < 0.0
true

```

That is, at `x` the function is changing sign.

From a mathematical perspective, a zero is guaranteed for a
*continuous* function. However, the computer algorithm doesn't assume
continuity, it just looks for changes of sign. As such, the algorithm
will  identify discontinuities, not just zeros. For example:

```jldoctest roots
julia> find_zero(x -> 1/x, (-1, 1))
0.0

```


The endpoints and function values can even be infinite:

```jldoctest roots
julia> find_zero(x -> Inf*sign(x), (-Inf, Inf))  # Float64 only
0.0

```


The basic algorithm used for bracketing when the values are simple
floating point values is a modification of the bisection method. For big float values, an
algorithm due to Alefeld, Potra, and Shi is used.

```jldoctest roots
julia> find_zero(sin, (big(3), big(4)))    # uses a different algorithm than for (3,4)
3.141592653589793238462643383279502884197169399375105820974944592307816406286198
```


The algorithms of Alefeld, Potra, and Shi and the well known algorithm of Brent, also start with a bracketing algorithm. For many problems these will take far fewer steps than the bisection algorithm to reach convergence. These may be called directly. For example,

```jldoctest roots
julia> find_zero(sin, (3,4), Roots.A42())
3.141592653589793
```

This takes ``9`` function evaluations, the default method takes ``53``. The method is specified above in the third positional argument by `Roots.A42()`. This method is not exported, so must be qualified.


By default, bisection will converge to machine tolerance. This may
provide more accuracy than desired. A tolerance may be specified to
terminate early, thereby utilizing fewer resources. For example, this
uses ``4`` steps to reach accuracy to $1/16$ (without specifying `xatol` it uses
``51`` steps):

```jldoctest roots
julia> rt = find_zero(sin, (3.0, 4.0), xatol=1/16)
3.125

julia> rt - pi
-0.016592653589793116

```


## Non-bracketing problems

Bracketing methods have guaranteed convergence, but in general may require
many more function calls than are otherwise needed to produce an answer.  If a
good initial guess is known, then the `find_zero` function provides an
interface to some different iterative algorithms that are more
efficient. Unlike bracketing methods, these algorithms may not
converge to the desired root if the initial guess is not well chosen.

The default algorithm is modeled after an algorithm used for
[HP-34 calculators](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf). This
algorithm is designed to be more forgiving of the quality of the
initial guess at the cost of possibly performing many more steps than
other algorithms, as if the algorithm encounters a bracket, a bracketing method
will be used.

For example, the answer to our initial problem is visibly seen from a
graph to be near 1. Given this,
the zero is found through:

```jldoctest roots
julia> f(x) = cos(x) - x
f (generic function with 1 method)


julia> x = find_zero(f , 1)
0.7390851332151607

julia> x, f(x)
(0.7390851332151607, 0.0)

```

For the polynomial $f(x) = x^3 - 2x - 5$, an initial guess of 2 seems reasonable:

```jldoctest roots
julia> f(x) = x^3 - 2x - 5
f (generic function with 1 method)

julia> x = find_zero(f, 2)
2.0945514815423265

julia> x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
(2.0945514815423265, -8.881784197001252e-16, -1.0)

```

For even more precision, `BigFloat` numbers can be used

```jldoctest roots
julia> x = find_zero(sin, big(3))
3.141592653589793238462643383279502884197169399375105820974944592307816406286198

julia> x, sin(x), x - pi
(3.141592653589793238462643383279502884197169399375105820974944592307816406286198, 1.096917440979352076742130626395698021050758236508687951179005716992142688513354e-77, 0.0)

```

### Higher order methods

The default call to `fzero` uses a first order method and then
possibly bracketing, which involves potentially many more function
calls. There may be times where a more efficient algorithm is
sought.  For such, a higher-order method might be better suited. There
are algorithms `Order1` (secant method), `Order2`
([Steffensen](http://en.wikipedia.org/wiki/Steffensen's_method)),
`Order5`, `Order8`, and `Order16`. The order 1 or 2 methods are
generally quite efficient. The even higher order ones are potentially
useful when more precision is used. These algorithms are accessed by
specifying the method after the initial starting point:

```jldoctest roots
julia> f(x) = 2x - exp(-x)
f (generic function with 1 method)

julia> x = find_zero(f, 1, Order1())      # also fzero(f, 1, order=1)
0.3517337112491958

julia> x, f(x)
(0.3517337112491958, -1.1102230246251565e-16)

```

The above makes $8$ function calls, to the $54$ made with `Order0`.

```jldoctest roots
julia> f(x) = (x + 3) * (x - 1)^2
f (generic function with 1 method)

julia> x = find_zero(f, -2, Order2())
-3.0

julia> x, f(x)
(-3.0, 0.0)

```


```jldoctest roots
julia> x = find_zero(f, 2, Order8())
1.0000000027152591

julia> x, f(x)
(1.0000000027152591, 2.949052856287529e-17)

```

The latter shows that zeros need not be simple zeros  to be found.
A simple zero, $c$,has $f(x) = (x-c) \cdot g(x)$ where $g(c) \neq 0$.
Generally speaking, non-simple zeros are
expected to take many more function calls, as the methods are no
longer super-linear. This is the case here, where `Order2` uses $51$
function calls, `Order8` uses $42$, and `Order0` takes  $80$. The `Roots.Order2B` method is useful
when a multiplicity is expected.

To investigate an algorithm and its convergence, the argument
`verbose=true` may be specified.


For some functions, adjusting the default tolerances may be necessary
to achieve convergence. The tolerances include `atol` and `rtol`, which are
used to check if $f(x_n) \approx 0$;
`xatol` and `xrtol`, to check if $x_n \approx x_{n-1}$; and `maxevals`  to limit the
number of steps in the algorithm.



## Classical methods

The package provides some classical methods for root finding:
`Roots.newton`, `Roots.halley`, and `Roots.secant_method`. (Currently
these are not exported, so must be prefixed with the package name to
be used.) We can see how each works on a problem studied by Newton
himself. Newton's method uses the function and its derivative:

```jldoctest roots
julia> f(x) = x^3 - 2x - 5
f (generic function with 1 method)

julia> fp(x) = 3x^2 - 2
fp (generic function with 1 method)

julia> x = Roots.newton(f, fp, 2)
2.0945514815423265

julia> x, f(x)
(2.0945514815423265, -8.881784197001252e-16)

```

To see the algorithm in progress, the argument `verbose=true` may be
specified.

Alternatively, `Roots.Newton()` can be specified as the method for `find_zero`. The
functions are specified using a tuple:

```jldoctest roots
julia> find_zero((f,fp), 2, Roots.Newton())
2.0945514815423265

```

The secant method typically needs two starting points, though a second
one is computed if only one is given. Here we start with 2 and 3,
specified through a tuple:

```jldoctest roots
julia> x = Roots.secant_method(f, (2,3))
2.094551481542327

julia> x, f(x)
(2.094551481542327, 3.552713678800501e-15)

```

Starting with a single point is also supported:

```jldoctest roots
julia> Roots.secant_method(f, 2)
2.0945514815423265

```

(This is like `Order1()`, or `Roots.Secant()`, but the implementation is
faster, as the framework is bypassed, and fewer checks on convergence
are used. This method can be used when speed is very important.)

Halley's method has cubic convergence, as compared to Newton's
quadratic convergence. It uses the second derivative as well:

```jldoctest roots
julia> fpp(x) = 6x
fpp (generic function with 1 method)

julia> x = Roots.halley(f, fp, fpp, 2)
2.0945514815423265

julia> x, f(x), sign(f(prevfloat(x)) * f(nextfloat(x)))
(2.0945514815423265, -8.881784197001252e-16, -1.0)

```

(Halley's method takes 3 steps, Newton's 4, but Newton's uses 5
function calls to Halley's 10.)

For many functions, their derivatives can be computed automatically. The
`ForwardDiff` package provides a means. Here we define an operator `D`
to compute a derivative:

```jldoctest roots
julia> using ForwardDiff

julia> function D(f, n::Int=1)
           n <= 0 && return f
           n == 1 && return x -> ForwardDiff.derivative(f,float(x))
           D(D(f,1),n-1)
       end
D (generic function with 2 methods)

julia> dfᵏs(f,k) = ntuple(i->D(f,i-1), Val(k+1)) # (f, f′, f′′, …)
dfᵏs (generic function with 1 method)
```


```jldoctest roots
julia> find_zero((f,D(f)), 2, Roots.Newton())
2.0945514815423265

```

Or, for Halley's method:

```jldoctest roots
julia> find_zero((f, D(f), D(f,2)), 2, Roots.Halley())
2.0945514815423265

```

The familiy of solvers implemented in `Roots.LithBoonkkampIJzerman(S,D)` where `S` is the number of prior points used to generate the next, and `D` is the number of derivatives used, has both the secant method (`S=2, D=0`) and Newton's method (`S=1, D=1`) as members, but also provides others. By adding more memory or adding more derivatives the convergence rate increases, at the expense of more complicated expressions or more function calls per step.

```
julia> find_zero(dfᵏs(f, 0), 2, Roots.LithBoonkkampIJzerman(3,0)) # like secant
2.0945514815423265

julia> find_zero(dfᵏs(f, 1), 2, Roots.LithBoonkkampIJzerman(2,1)) # like Newton
2.0945514815423265

julia> find_zero(dfᵏs(f, 2), 2, Roots.LithBoonkkampIJzerman(2,2)) # like Halley
2.0945514815423265
```

## The problem-algorithm-solve interface

The problem-algorithm-solve interface is a pattern popularized in `Julia` by the `DifferentialEquations.jl` suite of packages. The pattern consists of setting up a *problem* then *solving* the problem by specifying an *algorithm*. This is very similar to what is specified in the `find_zero(f, x0, M)` interface where `f` and `x0` specify the problem, `M` the algorithm, and `find_zero` calls the solver.

To break this up into steps, `Roots` has methods `ZeroProblem` and `init`, `solve`, and `solve!` from the `CommonSolve.jl` package.

Consider solving ``\sin(x) = 0`` using the `Secant` method starting with the interval ``[3,4]``.

```jldoctest roots
julia> f(x) = sin(x)
f (generic function with 1 method)

julia> x0 = (3, 4)
(3, 4)

julia> M = Roots.Secant()
Roots.Secant()

julia> Z = ZeroProblem(f, x0)
ZeroProblem{typeof(f), Tuple{Int64, Int64}}(f, (3, 4))

julia> solve(Z, M)
3.141592653589793
```

Changing the method is easy:

```jldoctest roots
julia> solve(Z, Roots.Order2())
3.1415926535897944
```

The `solve` interface works with parameterized functions, as well:

```jldoctest roots
julia> g(x,p) = cos(x) - x/p
g (generic function with 1 method)

julia> Z = ZeroProblem(g, (0.0, pi/2))
ZeroProblem{typeof(g), Tuple{Float64, Float64}}(g, (0.0, 1.5707963267948966))

julia> solve(Z, Roots.Secant(), 2) # p=2
1.0298665293222589

julia> solve(Z, Bisection(), 3)
1.1701209500026262
```
