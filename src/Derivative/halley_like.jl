struct HalleyState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δ::T
    ΔΔ::T
    fxn1::S
    fxn0::S
end

# we compute one step here to get x₁
function init_state(M::AbstractHalleyLikeMethod, F::Callable_Function, x)
    x₀ = float(first(x))
    T = eltype(x₀)
    fx₀, (Δ, ΔΔ) = F(x₀)
    Δx = calculateΔ(M, Δ, ΔΔ)
    x₁::T = x₀ - Δx
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

function init_state(::AbstractHalleyLikeMethod, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    fx₁, (Δ::T, ΔΔ::T) = F(x₁)
    HalleyState(x₁, x₀, Δ, ΔΔ, fx₁, fx₀)
end

initial_fncalls(M::AbstractHalleyLikeMethod) = 2 * 3
fn_argout(::AbstractHalleyLikeMethod) = 3

function update_state(
    M::AbstractΔMethod,
    F,
    o::HalleyState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1::T, r2::T = o.Δ, o.ΔΔ

    Δ::T = calculateΔ(M, r1, r2)
    if isissue(Δ)
        log_message(l, "Issue with computing `Δ`")
        return (o, true)
    end

    xn1::T = xn - Δ
    fxn1::S, (r1, r2) = F(xn1)
    incfn(l, 3)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1
    @set! o.Δ = r1
    @set! o.ΔΔ = r2

    return o, false
end

"""
    Roots.Halley()

Implements Halley's [method](https://en.wikipedia.org/wiki/Halley%27s_method), `xᵢ₊₁ = xᵢ
- (f/f')(xᵢ) * (1 - (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)^(-1)` This method
is cubically converging, it requires ``3`` function calls per
step. Halley's method finds `xₙ₊₁` as the zero of a hyperbola at the
point `(xₙ, f(xₙ))` matching the first and second derivatives of `f`.

The function, its derivative and second derivative can be passed either as a tuple of ``3`` functions *or*
as a function returning values for ``(f, f/f', f'/f'')``, which could be useful when function evaluations are expensive.

## Examples

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley()) ≈ π
true

julia> function f(x)
       s,c = sincos(x)
       (s, s/c, -c/s)
       end;

julia> find_zero(f, 3.0, Roots.Halley()) ≈ π
true
```

This can be advantageous if the derivatives are easily computed from
the computation for f, but otherwise would be expensive to compute separately.

----

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ ≈ -(2f'⋅f''' -3⋅(f'')²)/(12⋅(f'')²) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct Halley <: AbstractΔMethod end

"""
    Roots.QuadraticInverse()

Implements the [quadratic inverse method](https://doi.org/10.2307/2322644) also known as
[Chebyshev's method](https://dl.acm.org/doi/10.1080/00207160802208358),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 + (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)`.
This method is cubically converging, it requires ``3`` function calls per step.

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

"""
CHEBYSHEV-LIKE METHODS AND QUADRATIC EQUATIONS (J. A. EZQUERRO, J. M. GUTIÉRREZ, M. A. HERNÁNDEZ and M. A. SALANOVA)
"""
struct ChebyshevLike <: AbstractΔMethod end

"""
An acceleration of Newton's method: Super-Halley method (J.M. Gutierrez, M.A. Hernandez
"""
struct SuperHalley <: AbstractΔMethod end

# also Euler. Fits a parabola to point (x_n, f(x_n))
struct IrrationalHalley <: AbstractΔMethod end

"""
    Roots.Schroder()


Schröder's method, like Halley's method, utilizes `f`, `f'`, and
`f''`. Unlike Halley it is quadratically converging, but this is
independent of the multiplicity of the zero (cf. Schröder, E. "Über
unendlich viele Algorithmen zur Auflösung der Gleichungen."
Math. Ann. 2, 317-365, 1870;
[mathworld](http://mathworld.wolfram.com/SchroedersMethod.html)).

Schröder's method applies Newton's method to `f/f'`, a function with
all simple zeros.


## Example

```
m = 2
f(x) = (cos(x)-x)^m
fp(x) = (-x + cos(x))*(-2*sin(x) - 2)
fpp(x) = 2*((x - cos(x))*cos(x) + (sin(x) + 1)^2)
find_zero((f, fp, fpp), pi/4, Roots.Halley())     # 14 steps
find_zero((f, fp, fpp), 1.0, Roots.Schroder())    # 3 steps
```

(Whereas, when `m=1`, Halley is 2 steps to Schröder's 3.)

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```
find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Schroder())
```

This can be advantageous if the derivatives are easily computed from
the value of `f`, but otherwise would be expensive to compute.

The error, `eᵢ = xᵢ - α`, is the same as `Newton` with `f` replaced by `f/f'`.

"""
struct Schroder <: AbstractΔMethod end
const Schroeder = Schroder # either spelling
const Schröder = Schroder

## r1, r2 are o.Δ, o.ΔΔ
calculateΔ(::IrrationalHalley, r1, r2) = 2 / (1 + sqrt(1 - 2r1 / r2)) * r1
calculateΔ(::Halley, r1, r2) = 2r2 / (2r2 - r1) * r1
calculateΔ(::QuadraticInverse, r1, r2) = (1 + r1 / (2r2)) * r1
calculateΔ(::ChebyshevLike, r1, r2) = (1 + r1 / (2r2) * (1 + r1 / r2)) * r1
calculateΔ(::SuperHalley, r1, r2) = (1 + r1 / (2r2 - 2r1)) * r1
calculateΔ(::Schroder, r1, r2) = r2 / (r2 - r1) * r1

## Bracketed versions of Halley and Chebyshev using Alefeld, Potra, Shi approach
## calculateΔ has *different* calling pattern
## May take more steps and function evaluations, but should always converge

"""
    BracketedHalley

For a bracket `[a,b]`, uses the [`Roots.Halley`](@ref) method starting at the `x` value for which `fa` or `fb` is closest to `0`. Uses the `Roots.AbstractAlefeldPotraShi` framework to enforce the bracketing, taking an additional double secant step each time.
"""
struct BracketedHalley <: AbstractAlefeldPotraShi end
fn_argout(::BracketedHalley) = 3
fncalls_per_step(::BracketedHalley) = 3
function calculateΔ(::BracketedHalley, F::Callable_Function, c, ps)
    a, b = ps.a, ps.b

    fc, (Δ, Δ₂) = F(c)
    d = calculateΔ(Halley(), Δ, Δ₂)
    !(a <= c - d <= b) && (d = Δ)        # Newton
    d, ps
end

"""
    BracketedChebyshev

For a bracket `[a,b]`, uses the [`Roots.QuadraticInverse`](@ref) method starting at the `x` value for which `fa` or `fb` is closest to `0`. Uses the `Roots.AbstractAlefeldPotraShi` framework to enforce the bracketing, taking an additional double secant step each time.
"""
struct BracketedChebyshev <: AbstractAlefeldPotraShi end
fn_argout(::BracketedChebyshev) = 3
fncalls_per_step(::BracketedChebyshev) = 3
function calculateΔ(::BracketedChebyshev, F::Callable_Function, c, ps)
    a, b = ps.a, ps.a
    fc, (Δ, Δ₂) = F(c)
    d = calculateΔ(QuadraticInverse(), Δ, Δ₂)
    !(a <= c - d <= b) && (d = Δ)        # Newton
    d, ps
end

"""
    BracketedSchroder

For a bracket `[a,b]`, uses the [`Roots.Schroder`](@ref) method starting at the `x` value for which `fa` or `fb` is closest to `0`. Uses the `Roots.AbstractAlefeldPotraShi` framework to enforce the bracketing, taking an additional double secant step each time.
"""
struct BracketedSchroder <: AbstractAlefeldPotraShi end
fn_argout(::BracketedSchroder) = 3
fncalls_per_step(::BracketedSchroder) = 3
function calculateΔ(::BracketedSchroder, F::Callable_Function, c, ps)
    a, b = ps.a, ps.b

    fc, (Δ, Δ₂) = F(c)
    d = calculateΔ(Schroder(), Δ, Δ₂)
    !(a <= c - d <= b) && (d = Δ)        # Newton
    d, ps
end
