# Root finding functions for Julia

This package contains simple routines for finding the roots of continuous scalar functions of a single real variable. Two types of routine are available: bracketing and derivative based. The bracketing algorithms are more robust and pretty efficient, however you need to supply points either side of the zero (the bracket). The derivative based methods have faster convergence and high accuracy, but are not guaranteed to converge. Additionally, you must supply a derivative (although use of the Calculus package can help with that).

## Usage examples

### Finding roots without derivatives

The function `fzero` should be your first choice is you know a bracket for the
zero:

```julia
using Roots

root = fzero(sin, -0.5, 0.5)
println(sin(root))

f(x) = (x - 1)^3
root = fzero(f, 0.5, 2.0)
println(f(root))
```

You can pass in an optional `tolerance` parameter to specify the desired
accuracy (the default is 0.0 which specifies machine precision), and a
`max_iter` parameter to limit the number of iterations:

```julia
using Roots

f(x) = 2x*exp(-20) - 2*exp(-20x) + 1

root = fzero(f, 0.0, 1.0; tolerance=1e-10, max_iter=100)
println(f(root))
```

### Finding roots with derivatives

If you know the derivative of your function, you can use the Newton method
with a single initial guess (if you don't know the derivative, then you can
try using the Calculus package, as described later):

```julia
using Roots

f(x) = exp(x) - cos(x)
fp(x) = exp(x) + sin(x)

root = newton(f, fp, 3.0)
println(f(root))
```

If you are lucky enough to know the second derivative, you can also use
Halley's method, which has cubic convergence (rather than quadratic for Newton):

```julia
using Roots

f(x) = exp(x) - cos(x)
fp(x) = exp(x) + sin(x)
fpp(x) = exp(x) + cos(x)

root = halley(f, fp, fpp, 3.0)
println(f(root))
```

### Finding roots with numerical derivatives

It is also possible to use the marvelous [Calculus package](https://github.com/johnmyleswhite/Calculus.jl) to estimate the derivative using finite differences:

```julia
using Roots
using Calculus

f(x) = exp(x) - cos(x)
root = newton(f, derivative(f), 3.0)
println(f(root))
```

The precision is reduced slightly in this case, so unless you have a good
reason not to use it, the `fzero` method should be preferred.

