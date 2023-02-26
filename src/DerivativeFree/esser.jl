### Order2B() Esser method
"""
    Roots.Order2B()
    Roots.Esser()

Esser's method. This is a quadratically convergent method that, like
Schroder's method, does not depend on the multiplicity of the
zero. Schroder's method has update step `x - r2/(r2-r1) * r1`, where `ri =
fⁱ⁻¹/fⁱ`. Esser approximates `f' ~ f[x-h, x+h], f'' ~
f[x-h,x,x+h]`, where `h = fx`, as with Steffensen's method, Requiring 3
function calls per step. The implementation `Order2B` uses a secant
step when `|fx|` is considered too large.


Esser, H. Computing (1975) 14: 367. DOI: [10.1007/BF02253547](https://doi.org/10.1007/BF02253547)
Eine stets quadratisch konvergente Modifikation des Steffensen-Verfahrens


## Examples

```
f(x) = cos(x) - x
g(x) = f(x)^2
x0 = pi/4
find_zero(f, x0, Order2(), verbose=true)        #  3 steps / 7 function calls
find_zero(f, x0, Roots.Order2B(), verbose=true) #  4 / 9
find_zero(g, x0, Order2(), verbose=true)        #  22 / 45
find_zero(g, x0, Roots.Order2B(), verbose=true) #  4 / 10
```
"""
struct Esser <: AbstractSecantMethod end
struct Order2B <: AbstractSecantMethod end

function update_state(
    M::Order2B,
    fs,
    o::AbstractUnivariateZeroState,
    options,
    l=NullTracks(),
)
    update_state_guarded(M, Secant(), Esser(), fs, o, options, l)
end

function update_state(
    ::Esser,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    x1, fx1 = o.xn1, o.fxn1

    f0 = fx1

    f1  = F(x1 + f0 * oneunit(T) / oneunit(S))
    f_1 = F(x1 - f0 * oneunit(T) / oneunit(S))
    incfn(l, 2)

    # h = f0
    # r1 = f/f' ~ f/f[x+h,x-h]
    # r2 = f'/f'' ~ f[x+h, x-h]/f[x-h,x,x+h]
    r1 = f0 * 2 * f0 / (f1 - f_1) * oneunit(T) / oneunit(S)
    r2 = (f1 - f_1) / (f1 - 2 * f0 + f_1) * f0 / 2 * oneunit(T) / oneunit(S)

    k = r2 / (r2 - r1)  # ~ m

    if abs(k) <= 1e-2 * oneunit(k)
        log_message(l, "Estimate for multiplicity had issues. ")
        return o, true
    end

    delta = k * r1

    if isissue(delta)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0::T, x1::T = x1, x1 - delta
    fx0::S, fx1::S = fx1, F(x1)
    incfn(l)

    o = _set(o, (x1, fx1), (x0, fx0))

    return o, false
end
