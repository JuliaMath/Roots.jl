## Thukral 3,4,5 (2) is Schroder
"""
    AbstractThukralBMethod

Abstract type for `ThukralXB` methods for `X` being `2`,`3`,`4`, or `5`.

These are a family of methods which are
* efficient (order `X`) for non-simple roots (e.g. `Thukral2B` is the `Schroder` method)
* take `X+1` function calls per step
* require `X` derivatives. These can be passed as a tuple of functions, `(f, f', f'', …)`, *or* as
a function returning the ratios: `x -> (f(x), f(x)/f'(x), f'(x)/f''(x), …)`.

## Examples

```julia
using ForwardDiff
Base.adjoint(f::Function)  = x  -> ForwardDiff.derivative(f, float(x))
f(x) = (exp(x) + x - 2)^6
x0 = 1/4
find_zero((f, f', f''), x0, Roots.Halley())               # 14 iterations; ≈ 48 function evaluations
find_zero((f, f', f''), big(x0), Roots.Thukral2B())       #  3 iterations; ≈ 9 function evaluations
find_zero((f, f', f'', f'''), big(x0), Roots.Thukral3B()) #  2 iterations; ≈ 8 function evaluations
```


## Reference

*Introduction to a family of Thukral ``k``-order method for finding multiple zeros of nonlinear equations*,
R. Thukral, JOURNAL OF ADVANCES IN MATHEMATICS 13(3):7230-7237, DOI: [10.24297/jam.v13i3.6146](https://doi.org/10.24297/jam.v13i3.6146).
"""
abstract type AbstractThukralBMethod <: AbstractHalleyLikeMethod end
initial_fncalls(M::AbstractThukralBMethod) = fn_argout(M)

struct Thukral2B <: AbstractThukralBMethod end
fn_argout(::Thukral2B) = 3

struct Thukral3B <: AbstractThukralBMethod end
fn_argout(::Thukral3B) = 4

struct Thukral4B <: AbstractThukralBMethod end
fn_argout(::Thukral4B) = 5

struct Thukral5B <: AbstractThukralBMethod end
fn_argout(::Thukral5B) = 6

struct ThukralBState{N,T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δs::NTuple{N,T}
    fxn1::S
    fxn0::S
end

function init_state(M::AbstractThukralBMethod, F::Callable_Function, x)
    x₁ = float(first(x))
    fx₁, Δs = F(x₁)
    state = init_state(M, F, nan(x₁), x₁, nan(fx₁), fx₁; Δs=Δs)
end

function init_state(
    M::AbstractThukralBMethod,
    F,
    x₀::T,
    x₁::T,
    fx₀,
    fx₁;
    Δs=nothing,
) where {T}
    ThukralBState(x₁, x₀, NTuple{fn_argout(M) - 1,T}(Δs), fx₁, fx₀)
end

function update_state(
    M::AbstractThukralBMethod,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    x₀ = o.xn1

    Δ = compute_thukral_Δ(M, o)

    isissue(Δ) && return (o, true)
    x₁::T = x₀ - Δ
    fx₁::S, Δs = F(x₁)
    incfn(l, fn_argout(M))

    @set! o.xn0 = x₀
    @set! o.fxn0 = o.fxn1
    @set! o.Δs = NTuple{fn_argout(M) - 1,T}(Δs)
    @set! o.xn1 = x₁
    @set! o.fxn1 = fx₁

    return (o, false)
end

function compute_thukral_Δ(M::Thukral2B, o)
    r₁, r₂ = o.Δs
    t₁, t₂ = 1 / r₁, 1 / r₂
    Δ = one(o.xn1)
    Δ /= (t₁ - t₂)
    Δ
end

function compute_thukral_Δ(M::Thukral3B, o)
    r₁, r₂, r₃ = o.Δs
    t₁, t₂, t₃ = 1 / r₁, 1 / r₂, 1 / r₃
    Δ = (2t₁ - 2t₂)
    Δ /= (2t₁^2 - 3t₁ * t₂ + t₂ * t₃)
    Δ
end

function compute_thukral_Δ(M::Thukral4B, o)
    r₁, r₂, r₃, r₄ = o.Δs
    t₁, t₂, t₃, t₄ = 1 / r₁, 1 / r₂, 1 / r₃, 1 / r₄
    Δ = 6t₁^2 - 9t₁ * t₂ + 3t₂ * t₃
    Δ /= 6t₁^3 - 12 * t₁^2 * t₂ + 4t₁ * t₂ * t₃ - t₂ * t₃ * t₄ + 3 * t₁ * t₂^2
    Δ
end

function compute_thukral_Δ(M::Thukral5B, o)
    r₁, r₂, r₃, r₄, r₅ = o.Δs
    t₁, t₂, t₃, t₄, t₅ = 1 / r₁, 1 / r₂, 1 / r₃, 1 / r₄, 1 / r₅
    Δ = 24 * t₁^3 - 48t₁^2 * t₂ + 16 * t₁ * t₂ * t₃ - 4 * t₂ * t₃ * t₄ + 12t₁ * t₂^2
    Δ /=
        24t₁^4 - 60t₁^3 * t₂ + 20 * t₁^2 * t₂ * t₃ - 5 * t₁ * t₂ * t₃ * t₄ + 30t₁^2 * t₂^2 -
        10 * t₁ * t₂^2 * t₃ + t₂ * t₃ * t₄ * t₅
    Δ
end
