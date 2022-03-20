"""
    Roots.Schroder()


Schröder's method, like Halley's method, utilizes `f`, `f'`, and
`f''`. Unlike Halley it is quadratically converging, but this is
independent of the multiplicity of the zero (cf. Schröder, E. "Über
unendlich viele Algorithmen zur Auflösung der Gleichungen."
Math. Ann. 2, 317-365, 1870;
[mathworld](http://mathworld.wolfram.com/SchroedersMethod.html)). Schröder's
method applies Newton's method to `f/f'`, a function with all
simple zeros.


Example
```
m = 2
f(x) = (cos(x)-x)^m
fp(x) = (-x + cos(x))*(-2*sin(x) - 2)
fpp(x) = 2*((x - cos(x))*cos(x) + (sin(x) + 1)^2)
find_zero((f, fp, fpp), pi/4, Roots.Halley())    # 14 steps
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
struct Schroder <: AbstractHalleyLikeMethod end
const Schroeder = Schroder # either spelling
const Schröder = Schroder

## Shroder-like methods
function init_state(M::Schroder, F::Callable_Function, x)
    x₀ = float(first(x))
    fx₀, (Δ, ΔΔ) = F(x₀)
    x₁ = x₀ - ΔΔ / (ΔΔ - Δ) * Δ # m*r1
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

function update_state(
    method::Schroder,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options::UnivariateZeroOptions,
    l=NullTracks(),
) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ = r2 / (r2 - r1) * r1  # m * r1

    if isissue(Δ)
        log_message(l, "Issue with increment")
        return o, true
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
