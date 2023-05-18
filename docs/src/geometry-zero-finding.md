# Visualizing a step for various zero-finding algorithms

We illustrate the geometry behind a single step of several different, non-bracketing, zero-finding algorithms, beginning with, perhaps, the most famous, Newton's method.

## Newton's method

We load the `Plots`, `ForwardDiff`, and `Roots` packages:

```@example geometry
using Plots, ForwardDiff,  Roots
Base.adjoint(f::Function)  = x  -> ForwardDiff.derivative(f, float(x)) # f' will compute derivative
```

A zero-finding algorithm solves ``f(x) = 0`` or possibly ``f(x,p) = 0`` for a value of ``x``. Here we discuss iterative algorithms which take one *or more* past steps to produce the next step. (That is ``x_{n+1} = F(x_n, x_{n-1}, ..., x_1, x_0)``, for some ``F`` representing the algorithm).

[Newton's Method](https://en.wikipedia.org/wiki/Newton%27s_method) is a zero-finding *iterative algorithm* easily introduced in an introductory calculus class once the concept of a *tangent line* is presented.

The value ``x_{n+1}`` is described as the *intersection point* of the ``x``-axis with the tangent line through ``(x_n, f(x_n))``. To be explicit, we substitute ``(x_{n+1},0)`` into the tangent line equation ``y = f(x_n) + f'(x_n)\cdot(x-x_n)``:

```math
0 = f(x_n) + f'(x_n) \cdot (x_{n+1} - x_n).
```

Solving gives the update formula:

```math
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}.
```

The convergence is not guaranteed for all initial guesses, ``x_0``, but for a *simple zero* of a continuously differentiable function ``f(x)`` there is **some** interval about the zero, ``\alpha``, such that *quadratic convergence* will happen.

The geometry of Newton's method can be illustrated by graphing the tangent line.

The function ``f(x) = x^5 - x - 1`` does not have a readily available closed-form solution for its lone real zero, it being a fifth-degree polynomial. However, a graph, or other means, can show the function has one zero between ``1`` and ``2``, closer to ``1``. Starting with ``x_0=1.4``, we get a visual of ``x_1`` as follows:

```@example geometry
f(x) = x^5 - x - 1
x0 = 1.4
α = find_zero((f, f'), x0, Roots.Newton())

tl(x) = f(x0) + f'(x0)*(x-x0)
x1 = x0 - f(x0)/f'(x0)

p = plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)

plot!(tl; color="red", linewidth=3)

scatter!([x0, x1], [0, 0]; markercolor="blue")
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom)])

scatter!([x0], [f(x0)]; markercolor=:blue)

scatter!([α], [0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```


We used `Roots.Newton()` to identify the zero.

## Secant method

The secant method is much older than Newton's method, though similar in that the intersection of a line with the ``x``-axis is used as the next step in the algorithm. The slope of the secant line is (historically) easy to compute, unlike the slope of the tangent line which requires the notion of a derivative. The secant method begins with *two* initial points, ``x_0`` and ``x_1`` and uses the secant line instead of the tangent line. The secant line has slope ``(f(x_1)-f(x_0))/(x_1-x_0)``. This yields the algorithm:

```math
x_{n+1} = x_n - \left(\frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}\right)^{-1} \cdot f(x_n).
```

We can visualize the secant method easily enough. Suppose we start with ``x_0=1.4`` and ``x_1=1.3``:

```@example geometry
x0, x1 = 1.4, 1.3
x2 = x1 - (x1-x0)/(f(x1)-f(x0)) * f(x1)
sl(x) = f(x1) + (f(x1)-f(x0))/(x1-x0) * (x-x1)

p = plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)

plot!(sl, color=:red, linewidth=3)

scatter!([x0, x1, x2], [0,0,0]; markercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom), (x2,0,"x2", :bottom)])
scatter!([x0, x1], [f(x0), f(x1)]; markercolor=:blue)
scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```

The secant method is implemented in `Secant()`.

Steffensen's method (`Root.Steffensen()`) is related to the secant method, though the points are not ``x_n`` and ``x_{n-1}``, rather  ``x_n + f(x_n)`` and ``x_n``. As ``x_n`` gets close to ``\alpha``, ``f(x_n)`` gets close to ``0``, so this method converges at an asymptotic rate like Newton's method. (Though with a tradeoff, as the secant method needs only one new function evaluation per step, Steffensen's require two.)

### An inverse view

The secant line is a natural object as presented above, but can be viewed a bit differently. Consider the two points ``(f(x_0), x_0)`` and ``(f(x_1),x_1)``. Two non-equal points determine a line. In this case, we have inverted the ``x`` and ``y`` values, so we invert the coordinates of the line. To find ``x = my + b``, or some other form of the line involves solving two equations with two unknowns. Each equation comes by using the known point:

```math
\begin{align*}
x_0 &= m \cdot f(x_0) + b\\
x_1 &= m \cdot f(x_1) + b
\end{align*}
```

This *linear* set of equations can be solved, some `SymPy` code would look like:

```
using SymPy
@syms x0, y0, x1, y1, m ,b
u = solve([x0 ~ y0 * m + b, x1 ~ y1 * m + b], (m,b))
```

Producing

```
Dict{Any, Any} with 2 entries:
  b => (-x0*y1 + x1*y0)/(y0 - y1)
  m => (x0 - x1)/(y0 - y1)
```

The value of `m` is the reciprocal of the slope, as we have inverted the perspective. The value of `b` is where the inverse line intersects the ``x`` axis and is the same as the secant method algorithm:

```
sm = x1 - y1 * (x1-x0)/(y1-y0)
simplify(sm - u[b])
```

leading to:

```
0
```

## Inverse quadratic and cubic methods

Brent's method (`Roots.Brent()`) is a bracketing method which utilizes an inverse quadratic step to speed up convergence beyond the secant method. The inverse quadratic step uses the fact that three (non-collinear) points determine a quadratic polynomial. As above, this is done with the inverse of the points, ``(x_{n-2}, f(x_{n-2}))``, ``(x_{n-1}, f(x_{n-1}))``, and ``(x_{n}, f(x_{n}))``. Using the same method illustrated above, it can be shown that, with ``\Delta_{i,j} = f(x_i) - f(x_j)``:

```math
x_{n+1} = \frac{x_{n-2}f(x_{n-1})f(x_n)}{\Delta_{i-2, i-1}\Delta_{i-2,i}}
+ \frac{f(x_{n-2})x_{n-1}f(x_n)}{\Delta_{n-1, n-2}\Delta_{n-1,n}}
+ \frac{f(x_{n-2})f(x_{n-1})x_n}{\Delta_{n,n-2}\Delta_{n,n-1}}.
```


```
x0, x1, x2 = xs = 1.4, 1.3, 1.2
fx0, fx1, fx2 = ys = f.(xs)

@syms x[0:2], y[0:2], a, b, c
u = solve([xᵢ ~ a*yᵢ^2 + b * yᵢ + c for (xᵢ, yᵢ) ∈ zip(x, y)], (a, b, c))
x3 = u[c]
for (k, v) ∈ u
  for (xᵢ, yᵢ, x,y) ∈ zip(x, y, xs, ys)
    v = v(xᵢ => x, yᵢ => y)
  end
  u[k] = v
end
u[a], u[b], u[c]
```

Which returns

```
(-0.00930682152998560, 0.104752944517765, 1.17057129242798)
```

This last value, `c`, is also computed in the `(3,0)` method of the `LithBoonkkampIJzerman` algorithm, which implements this method:

```@example geometry
x0, x1, x2 = 1.4, 1.3, 1.2
xs = [x0, x1,x2]
x3 = Roots.lmm(Roots.LithBoonkkampIJzerman{3,0}(), xs, f.(xs))
```

With this, we can visualize:

```@example geometry
a, b, c = (-0.00930682152998560, 0.104752944517765, 1.17057129242798)
iq(y) = a * y^2 + b * y + c

p = plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)

ys′ = range(f(1.1), f(1.5), length=100)
plot!(iq.(ys′), ys′; color=:red, linewidth=3)

scatter!(xs, f.(xs);  markercolor=:blue)
scatter!([x3], [f(x3)]; markercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom),
	(x2,0,"x2", :bottom),  (x3,0,"x3", :bottom)])
scatter!(xs, zero.(xs);  markercolor=:blue)

scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```

Inverse cubic is similar to the above, though we have ``4`` past points used to determine the next one. In this example, we solve the resulting linear system of equations numerically (the `Roots.LithBoonkkampIJzerman{4,0}()` method implements the algorithm):

```@example geometry
xs = [1.4, 1.35, 1.3, 1.25]
ys = f.(xs)
A = zeros(Float64, 4, 4)
for i ∈ reverse(0:3)
    A[:,4-i] .= ys.^i
end
a, b, c, d = A \ xs
ic(y) = a * y^3 + b * y^2 + c * y + d

p = plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)
plot!(ic.(ys′), ys′;  color=:red, linewidth=3)

x4 = d
scatter!(vcat(xs, x4), zeros(5); markercolor=:blue)
for (i,x) ∈ enumerate(xs)
    annotate!([(x, 0, "x$(i-1)", :bottom)])
end
annotate!([(x4, 0, "x4", :bottom)])
scatter!(xs, f.(xs); markercolor=:blue)
scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```


We can see graphically that for this function and the chosen values, the inverse quadratic and inverse cubic fits are very close to the actual zero, suggesting a rapid convergence. `Roots.Brent()`, `Roots.Chandrapatlu()`, and `Roots.AlefeldPotraShi()` are different bracketing methods which use an inverse quadratic step *when* the algorithm deems it productive, falling back to other estimates when not. Similarly, the inverse cubic step is utilized by `Roots.A42()`, as possible. The `LithBoonkkampIJzerman{S,0}` methods use ``S`` previous points (``S \geq 2``) and the corresponding inverse polynomial step to progress. Since these are not bracketed, the algorithms are only guaranteed to converge for nearby initial guesses.


## Higher derivative variations on Newton's method

We can visualize Newton's method differently than as an intersection of the ``x``-axis with a specific tangent line, rather we can think of it as an intersection of the ``x``-axis with *a* line ``ax + by = c`` with *two* constraints:

* the point ``(x_n, f(x_n))`` is on the line
* the line matches the slope of the tangent line at that point (a tangency condition)

Combined, these say that the constrained line has slope ``f'(x_n)`` and goes through the point ``(x_n, f(x_n))``, so is the tangent line; Newton's method follows.

When a quadratic approximation to the graph of ``f(x)`` is chosen *at* ``(x_n, f(x_n))`` other algorithms become possible.

[Geometric constructions of iterative functions to solve nonlinear equations](https://doi.org/10.1016/S0377-0427(03)00420-5) by Amat, Busquier, and Gutiérrez has a systematic approach we follow.

### Euler method

Consider the quadratic expression ``y + ax^2 + bx + c = 0``.
Assuming the expression goes through the point ``(x_n, f(x_n))`` and the tangency conditions: ``y'(x) = f'(x)`` **and** ``y''(x) = f''(x)``, we get the second-order Taylor polynomial ``y(x) = f(x_n) + f'(x_n)(x-x_n) + f''(x_n)/2 \cdot(x-x_n)^2`` as the solution.

Let

```math
L_f(x) = \frac{f(x)/f'(x)}{f'(x)/f''(x)} = \frac{f(x)f''(x)}{f'(x)^2}.
```

Then the intersection of the curve with the ``x``-axis can be represented as:

```math
x_{n+1} = x_n - \frac{2}{1 + \sqrt{1 - 2L_f(x_n)}} \cdot \frac{f(x_n)}{f'(x_n)}.
```

This is known as Euler's method or the irrational Halley method and implemented in `Roots.IrrationalHalley()`:

```@example geometry
L_f(x) = f(x) * f''(x) / (f'(x))^2

x0 = 1.4
x1 = x0 - 2 / (1 + sqrt(1 - 2L_f(x0))) * f(x0)/f'(x0)
t2(x) = f(x0) + f'(x0)*(x-x0) + f''(x0)/2 * (x - x0)^2


p = plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)
plot!(t2; color=:red, linewidth=3)

scatter!([x0, x1], [0,0]; markercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom)])
scatter!([x0], [f(x0)]; markercolor=:blue)
scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```


### Halley's method

The general form of a quadratic equation ``Ax^2 + By^2 + Cxy + Dx + Ey + F = 0`` is specialized above by setting ``B=C=0`` and ``E=1`` and then imposing the point ``(x_n, f(x_n))`` as a solution along with tangency conditions. The famous Halley's method can be seen as the specialization to a hyperbola: ``axy + y + bx + c = 0``. This yields the curve

```math
y  - (f(x_n)  + f'(x_n)(x-x_n) + \frac{f''(x_n)}{2f'(x_n)}(x-x_n) \cdot (y-f(x_n))) = 0
```

and the iterative algorithm:

```math
x_{n+1} = x_n - \frac{2}{2 - L_f(x_n)} \cdot \frac{f(x_n)}{f'(x_n)}.
```

We can visualize, as follows, using a contour plot to represent the hyperbola.

```@example geometry
x1 = x0 - 2 / (2 - L_f(x0)) * f(x0)/f'(x0)
F(x,y) = y - f(x0) - f'(x0)*(x-x0) - f''(x0)/(2f'(x0)) * (x-x0) * (y-f(x0))


plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)
xs′, ys′ = range(1.1, 1.5, length=50), range(f(1.1), f(1.5), length=50);
zs = [F(x,y) for y ∈ ys′, x ∈ xs′];
contour!(xs, ys, zs; levels = [0], color=:red, linewidth=3)

scatter!([x0, x1], [0,0]; markercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom)])

scatter!([x0], [f(x0)]; markercolor=:blue)

scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```

`Roots.Halley()` provides an implementation.

### Chebyshev's method

Chebyshev's method uses an inverse quadratic fit, specialized with ``ay^2 + y + bx + c = 0``, to compute the next iterate. It can be expressed via:

```math
\frac{-f''(x_n)}{2f'(x_n)^2}(y - f(x_n))^2 + y - f(x_n) - f'(x_n)(x-x_n) = 0
```

and the algorithm becomes:

```math
x_{n+1} = x_n - (1 + \frac{1}{2} L_f(x_n)) \frac{f(x_n)}{f'(x_n)}.
```

This is visualized in a similar manner as the last example:

```@example geometry
x1 = x0 + (1 - 1/2 * L_f(x0)) * f(x0) / f'(x0)
F(x, y) = -f''(x0)/(2f'(x0)^2) * (y-f(x0))^2 + y - f(x0)  - f'(x0) * (x- x0)

plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)
xs′, ys′ = range(1.1, 1.5, length=50), range(f(1.1), f(1.5), length=50);
zs = [F(x,y) for y ∈ ys′, x ∈ xs′];
contour!(xs, ys, zs; levels = [0], color=:red, linewidth=3)

scatter!([x0, x1], [0,0]; markermarkercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom)])

scatter!([x0], [f(x0)]; markermarkercolor=:blue)

scatter!([α],[0]; markermarkercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```

`Roots.InverseQuadratic()` provides an implementation; `Roots.ChebyshevLike()` an accelerated version.

Amat, Busquier, and Gutierrez also consider the hyperbola

```math
ay^2 + bxy + y + cx + d = 0
```

As there are ``4`` unknowns and only ``3`` constraints, the solution will depend on a parameter, they call ``b_n,`` yielding:

```math
\begin{align*}
c_n &= -f'(x_n)\\
a_n &= -\frac{f''(x_n)}{2f'(x_n)^2} - \frac{b_n}{f'(x_n)}\\
0 &= x - x_n  + (y-f(x_n)) \frac{1 + a_n \cdot (y - f(x_n))}{b_n\cdot (y-f(x_n)) + c_n}
\end{align*}
```

This gives the algorithm:

```math
x_{n+1} = x_n - \left(1 + \frac{1}{2}\frac{L_f(x_n)}{1 + b_n \cdot (f(x_n)/f'(x_n))}\right) \cdot \frac{f(x_n)}{f'(x_n)}.
```

Newton's method is recovered by letting ``b_n \rightarrow 0``, Chebyshev's method is when ``b_n=0``, Halley's method is when ``b_n = -f''(x_n)/(2f'(x_n))``.

The super-Halley method is when ``b_n = -f''(x_n)/f'(x_n)``. We can visualize this:

```@example geometry
cn = -f(x0)
bn = -f''(x0)/f'(x0)
an = -f''(x0)/(2f'(x0)) - bn/f'(x0)
x1 = x0 - (1 + 1/2 * (L_f(x0) / (1 + bn * (f(x0)/f'(x0))))) * f(x0)/f'(x0)
F(x, y) = x - x0 + (y-f(x0)) * (1 + an * (y - f(x0))) / (bn * (y - f(x0)) + cn)

plot(f, 1.1, 1.5; legend=false, linewidth=3)
plot!(zero)
xs′, ys′ = range(1.1, 1.5, length=50), range(f(1.1), f(1.5), length=50);
zs = [F(x,y) for y ∈ ys′, x ∈ xs′];
contour!(xs, ys, zs; levels = [0], color=:red, linewidth=3)

scatter!([x0, x1], [0,0]; markercolor=:blue)
annotate!([(x0,0,"x0", :bottom), (x1, 0, "x1", :bottom)])
scatter!([x0], [f(x0)]; markercolor=:blue)
scatter!([α],[0]; markercolor=:blue)
annotate!([(α, 0, "α", :top)])
p
```

`Roots.SuperHalley()` provides an implementation.

The author's discuss using different osculating curves, such as a cubic equation. As they mention, all their methods take the form (using [big O notation](https://en.wikipedia.org/wiki/Big_O_notation)):

```math
x_{n+1} = x_n - (1 + \frac{1}{2}L_f(x_n) + \mathcal{O}(L_f(x_n)^2)) \cdot \frac{f(x_n)}{f'(x_n)}.
```

The algorithms implemented in the `Roots.LithBoonkkampIJzerman{S,D}()` methods use a differential equations approach to evaluate the inverse of $f(x)$ at $0$. The methods all find inverse polynomial approximations to $f^{-1}$. The methods for $D=0$, are, as seen, for $S=2$ an inverse secant line, for $S=3$ an inverse quadratic approximation, for $S=4$ an inverse cubic; when $S=1$ and $D=1$ Euler's method turns into Newton's method, for $S=1$ and $D=2$ the inverse quadratic or Chebyshev method is used.
