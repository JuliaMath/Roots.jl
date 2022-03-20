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
    x₁::T = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

function init_state(::AbstractHalleyLikeMethod, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    fx₁, (Δ::T, ΔΔ::T) = F(x₁)
    HalleyState(x₁, x₀, Δ, ΔΔ, fx₁, fx₀)
end

initial_fncalls(M::AbstractHalleyLikeMethod) = 2 * 3
fn_argout(::AbstractHalleyLikeMethod) = 3


function update_state(
    method::AbstractΔMethod,
    F,
    o::HalleyState{T,S},
    options::UnivariateZeroOptions,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ::T = calculateΔ(method, r1, r2)
    if isissue(Δ)
        log_message(l, "Issue with computing `Δ`")
        return (o, true)
    end

    xn1::T = xn - Δ
    fxn1::S, (r1::T, r2::T) = F(xn1)
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

Implements Halley's [method](http://tinyurl.com/yd83eytb),
`xᵢ₊₁ = xᵢ - (f/f')(xᵢ) * (1 - (f/f')(xᵢ) * (f''/f')(xᵢ) * 1/2)^(-1)`
This method is cubically converging, but requires more function calls per step (3) than
other methods.

Example

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin, cos, x->-sin(x)), 3.0, Roots.Halley()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which
returns `(f, f/f',f'/f'')` as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x), -cos(x)/sin(x)), 3.0, Roots.Halley()) ≈ π
true
```

This can be advantageous if the derivatives are easily computed from
the computation for f, but otherwise would be expensive to compute separately.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₁ ≈ -(2f'⋅f''' -3⋅(f'')²)/(12⋅(f'')²) ⋅ eᵢ³` (all evaluated at `α`).

"""
struct Halley <: AbstractΔMethod end


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


"""
CHEBYSHEV-LIKE METHODS AND QUADRATIC EQUATIONS (J. A. EZQUERRO, J. M. GUTIÉRREZ, M. A. HERNÁNDEZ and M. A. SALANOVA)
"""
struct ChebyshevLike <: AbstractΔMethod end

"""
An acceleration of Newton's method: Super-Halley method (J.M. Gutierrez, M.A. Hernandez
"""
struct SuperHalley <: AbstractΔMethod end

## r1, r2 are o.Δ, o.ΔΔ
calculateΔ(method::Halley, r1, r2) = 2r2 / (2r2 - r1) * r1
calculateΔ(method::QuadraticInverse, r1, r2) = (1 + r1 / (2r2)) * r1
calculateΔ(method::ChebyshevLike, r1, r2) = (1 + r1 / (2r2) * (1 + r1 / r2)) * r1
calculateΔ(method::SuperHalley, r1, r2) = (1 + r1 / (2r2 - 2r1)) * r1
