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

To illustrate, consider the function $f(x) = \cos(x) - x$. From the
graph we see readily that $[0,1]$ is a bracket (which we emphasize
with an overlay):

```@example roots
using Plots, Roots
f(x) = cos(x) - x
plot(f, -2, 2)
plot!([0,1], [0,0], linewidth=2)
savefig("plot-f-1.svg"); nothing #hide
```

![](plot-f-1.svg)
