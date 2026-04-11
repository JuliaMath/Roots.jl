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

struct ThukralB{N} <: AbstractThukralBMethod end
const Thukral2B = ThukralB{2}
const Thukral3B = ThukralB{3}
const Thukral4B = ThukralB{4}
const Thukral5B = ThukralB{5}

initial_fncalls(::ThukralB{N}) where {N} =  N
fn_argout(::ThukralB{N}) where N = N + 1

struct ThukralBState{N,T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δs::NTuple{N,T}
    fxn1::S
    fxn0::S
end

function init_state(M::ThukralB{N}, F::Callable_Function, x) where N
    x₁ = float(first(x))
    fx₁, Δs = F(x₁)
    S = eltype(fx₁)
    state = init_state(M, F, nan(x₁), x₁, nan(fx₁), fx₁; Δs = ntuple(i -> S(Δs[i]), Val(N)))
end

function init_state(
    M::ThukralB{N},
    F,
    x₀::T,
    x₁::T,
    fx₀::S,
    fx₁::S;
    Δs=NTuple{N,S}(),
) where {T,S,N}
    x1, x0 = promote(x₁, x₀)
    fx1, fx0 = promote(fx₁, fx₀)
    ThukralBState(x1, x0, Δs, fx1, fx0)
end

function update_state(
    M::ThukralB{N},
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S,N}
    x₀ = o.xn1

    Δ = compute_thukral_Δ(M, o)

    isissue(Δ) && return (o, true)
    x₁::T = x₀ - Δ
    fx₁::S, Δs = F(x₁)
    incfn(l, fn_argout(M))

    @reset o.xn0 = x₀
    @reset o.fxn0 = o.fxn1
    @reset o.Δs = ntuple(i -> T(Δs[i]), Val(N))
    @reset o.xn1 = x₁
    @reset o.fxn1 = fx₁

    return (o, false)
end

function compute_thukral_Δ(M::ThukralB{2}, o)
    r₁, r₂ = o.Δs
    t₁, t₂ = 1 / r₁, 1 / r₂
    Δ = one(o.xn1)
    Δ /= (t₁ - t₂)
    Δ
end

function compute_thukral_Δ(M::ThukralB{3}, o)
    r₁, r₂, r₃ = o.Δs
    t₁, t₂, t₃ = 1 / r₁, 1 / r₂, 1 / r₃
    Δ = (2t₁ - 2t₂)
    Δ /= (2t₁^2 - 3t₁ * t₂ + t₂ * t₃)
    Δ
end

function compute_thukral_Δ(M::ThukralB{4}, o)
    r₁, r₂, r₃, r₄ = o.Δs
    t₁, t₂, t₃, t₄ = 1 / r₁, 1 / r₂, 1 / r₃, 1 / r₄
    Δ = 6t₁^2 - 9t₁ * t₂ + 3t₂ * t₃
    Δ /= 6t₁^3 - 12 * t₁^2 * t₂ + 4t₁ * t₂ * t₃ - t₂ * t₃ * t₄ + 3 * t₁ * t₂^2
    Δ
end

function compute_thukral_Δ(M::ThukralB{5}, o)
    r₁, r₂, r₃, r₄, r₅ = o.Δs
    t₁, t₂, t₃, t₄, t₅ = 1 / r₁, 1 / r₂, 1 / r₃, 1 / r₄, 1 / r₅
    Δ = 24 * t₁^3 - 48t₁^2 * t₂ + 16 * t₁ * t₂ * t₃ - 4 * t₂ * t₃ * t₄ + 12t₁ * t₂^2
    Δ /=
        24t₁^4 - 60t₁^3 * t₂ + 20 * t₁^2 * t₂ * t₃ - 5 * t₁ * t₂ * t₃ * t₄ + 30t₁^2 * t₂^2 -
        10 * t₁ * t₂^2 * t₃ + t₂ * t₃ * t₄ * t₅
    Δ
end
