struct ScalingFactorRegulaFalsi{R} <: AbstractRegulaFalsiMethod end

"""

    ScalingFactorRegulaFalsi(M::Symbol=AndersonBjork)
    RegulaFalsi{M}
    SCRF{M}

Implements several different *scaling factor* regular falsi methods following Diaz and Perez. The regular falsi method uses the secant method to find the next value with a bracketing interval providing the history. This is a linearly convergent method. By scaling, or adjusting, the ``f(x_i)`` values, superlinear convergence can be achieved.

The scaling factor refers to the ``f(x_i)`` values being scaled by `γ` which is dynamically computed based on ``\\xi=f(x_{i+1})/f(x_i)`` and ``\\zeta = -f(x_{i+1}) / f(x_{i-1})``.

The value for `M` is in `(:classic, :Illinois, :GIll01, :Pegasus, :AndersonBjork, :AB_GIll01, :Ford3, :Ford4,)`. The defaults is :AndersonBjork.



Examples
```
find_zero(x -> x^5 - x - 1, (-2, 2), ScalingFactorRegulaFalsi())
find_zero(x -> x^5 - x - 1, (-2, 2), ScalingFactorRegulaFalsi(:Illinois))
```

New scaling factors can be introduced by defining a method for `Fᵧ(::ScalingFactorRegulaFalsi{M},  ξ::T, ζ::T) where T` and then calling with `ScalingFactorRegulaFalsi{M}` for `M` a symbol.

Convergence rates from the reference paper are:

1. `:classic`---Classic Regula Falsi: p = 1 (linear).
2. `:Illinois`---Illinois method: p = 1.442.
3. `:Pegasus`---Pegasus method: p = 1.642.
4. `:AndersonBjork`---A&B method: 1.681 < p < 1.710.
5. `:Ford4`---Ford fourth method: p = 1.681.


References:

*A common framework for modified Regula Falsi methods and new methods of this kind* by Julio M. Fernández-Díaz, César O. Menéndez-Pérez [url](https://digibuo.uniovi.es/dspace/bitstream/handle/10651/66581/1-s2.0-S0378475422004335-main.pdf)

Note: Compared to similar methods in `FalsePosition` the parameterization used here has fewer issues with floating point differencees.

"""
ScalingFactorRegulaFalsi
ScalingFactorRegulaFalsi(x=:AndersonBjork) = ScalingFactorRegulaFalsi{x}()
const RegulaFalsi{R} = ScalingFactorRegulaFalsi{R}
const SCRF{R} = ScalingFactorRegulaFalsi{R}

# for testing, this might be helpful
_regula_falsi_names = (:classic, :Illinois, :GIll01, :Pegasus, :AndersonBjork, :AB_GIll01, :Ford3, :Ford4)

# factors
Fᵧ(::ScalingFactorRegulaFalsi{:classic}, ξ::T, ζ::T) where T = one(T)
Fᵧ(::ScalingFactorRegulaFalsi{:Illinois}, ξ::T, ζ::T) where T = one(T)/2
Fᵧ(::ScalingFactorRegulaFalsi{:GIll01}, ξ::T, ζ::T) where T = one(T)/10 # Generalized Illinois method with γ = 0.1
Fᵧ(::ScalingFactorRegulaFalsi{:Pegasus}, ξ::T, ζ::T) where T = one(T) / (one(T) + ξ)
Fᵧ(::ScalingFactorRegulaFalsi{:AndersonBjork}, ξ::T, ζ::T) where T = ξ < 1 ? (one(T) - ξ) : one(T)/2
Fᵧ(::ScalingFactorRegulaFalsi{:AB_GIll01}, ξ::T, ζ::T) where T = max(one(T) - ξ,one(T)/10) # AB + Generalized Illinois method with γ = 0.1
Fᵧ(::ScalingFactorRegulaFalsi{:Ford3},  ξ::T, ζ::T) where T = (one(T) - ξ + ζ) / (one(T) - ζ)
Fᵧ(::ScalingFactorRegulaFalsi{:Ford4},  ξ::T, ζ::T) where T = one(T) - ξ + ζ

function init_state(M::AbstractRegulaFalsiMethod, F, x₀::S, x₁::S, fx₀::T, fx₁::T) where {S, T}
    assert_bracket(fx₀,fx₁)

    c::T = (x₀ * fx₁ - x₁ * fx₀) / (fx₁ - fx₀)
    fc::S = F(c)

    if sign(fx₀) * sign(fc) < 0
        a, fa = x₀, fx₀
    else
        a, fa = x₁, fx₁
    end
    UnivariateZeroState(c, a, fc, fa)
end
initial_fncalls(M::AbstractRegulaFalsiMethod) = 3

function default_tolerances(
    ::ScalingFactorRegulaFalsi{K},
    ::AbstractUnivariateZeroState{T,S},
) where {K, T,S}
    xatol = eps(T)^3 * oneunit(T)
    xrtol = eps(T) * one(T) # unitless
    atol = 0 * oneunit(S)
    rtol = 0 * one(S)
    maxiters = 250
    strict = false
    (xatol, xrtol, atol, rtol, maxiters, strict)
end


function update_state(
    M::ScalingFactorRegulaFalsi,
    fs,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}

    xₙ₋₁::T, xₙ::T = o.xn0, o.xn1
    f̂xₙ₋₁::S, f̂xₙ::S = o.fxn0, o.fxn1

    c::T = (xₙ₋₁ * f̂xₙ - xₙ * f̂xₙ₋₁) / (f̂xₙ - f̂xₙ₋₁)
    if (c == xₙ₋₁ || c == xₙ) # try midpoint o/w this is stuck
        c = xₙ₋₁/2 + xₙ/2
    end
    fc::S = fs(c)

    incfn(l)


    iszero(fc) && return (_set(o, (c, fc)), true)

    xₙ₊₁, f̂xₙ₊₁ = c, fc

    if sign(f̂xₙ) * sign(fc) < 0
        xₙ, f̂xₙ = xₙ, f̂xₙ
    elseif sign(f̂xₙ) * sign(fc) > 0
        ξ::S = fc / f̂xₙ     # \xi
        ζ::S = - fc / f̂xₙ₋₁ # \zeta
        γ::S =  Fᵧ(M, ξ, ζ)
        γ =  max(zero(ξ), min(one(ξ), γ))

        xₙ = xₙ₋₁
        f̂xₙ = γ * f̂xₙ₋₁
        iszero(f̂xₙ) && (f̂xₙ = f̂xₙ₋₁/2) # Some engineering to avoid early termination by f̂(xₙ) == 0.0
    end

    o = _set(o, (xₙ₊₁, f̂xₙ₊₁), (xₙ, f̂xₙ))

    return (o, false)
end
