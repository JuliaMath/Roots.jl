"""
    Roots.QuadraticInverse()

Implements the [quadratic inverse method](https://doi.org/10.2307/2322644) also known as [Chebyshev's method]((https://dl.acm.org/doi/10.1080/00207160802208358)),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 + (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)`.
This method is cubically converging, but requires more function calls per step (3) than
other methods.

Example

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.QuadraticInverse()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.QuadraticInverse()) ≈ π
true
```

This can be advantageous if the derivatives are easily computed from
the computation for f, but otherwise would be expensive to compute separately.

The error, `eᵢ = xᵢ - α`, [satisfies](https://dl.acm.org/doi/10.1080/00207160802208358)
`eᵢ₊₁ ≈ (1/2⋅(f''/f')² - 1/6⋅f'''/f')) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct QuadraticInverse <: AbstractΔMethod end

calculateΔ(method::QuadraticInverse, r1, r2) = (1 + r1 / (2r2)) * r1

"""
    Roots.quadratic_inverse(f, fp, fpp, x0; kwargs...)

Implementation of the quadratic inverse method (cf `?Roots.QuadraticInverse()`).

Arguments:

* `f::Function` -- function to find zero of

* `fp::Function` -- derivative of `f`.

* `fpp:Function` -- second derivative of `f`.

* `x0::Number` -- initial guess

With the `FowardDiff` package derivatives may be computed automatically. For example,  defining
`D(f) = x -> ForwardDiff.derivative(f, float(x))` allows `D(f)` and `D(D(f))` to be used for the first and second
derivatives, respectively.

Keyword arguments are passed to `find_zero` using the `Roots.QuadraticInverse()` method.


"""
quadratic_inverse(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, QuadraticInverse(); kwargs...)

"""
CHEBYSHEV-LIKE METHODS AND QUADRATIC EQUATIONS (J. A. EZQUERRO, J. M. GUTIÉRREZ, M. A. HERNÁNDEZ and M. A. SALANOVA)
"""
struct ChebyshevLike <: AbstractΔMethod end

calculateΔ(method::ChebyshevLike, r1, r2) = (1 + r1 / (2r2) * (1 + r1 / r2)) * r1

chebyshev_like(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, ChebyshevLike(); kwargs...)

"""
An acceleration of Newton's method: Super-Halley method (J.M. Gutierrez, M.A. Hernandez
"""
struct SuperHalley <: AbstractΔMethod end

calculateΔ(method::SuperHalley, r1, r2) = (1 + r1 / (2r2 - 2r1)) * r1

superhalley(f, fp, fpp, x0; kwargs...) =
    find_zero((f, fp, fpp), x0, SuperHalley(); kwargs...)
