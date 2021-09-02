######
# methods that converge quickly even when α is non-simple
# These usually take 1 additional function call compared to a similar ordered
# method
######


"""
    Roots.Order1B()
    Roots.King()

A superlinear (order `1.6...`) modification of the secant method for multiple roots.
Presented in A SECANT METHOD FOR MULTIPLE ROOTS, by RICHARD F. KING, BIT 17 (1977), 321-328

The basic idea is similar to Schroder's method: apply the secant method
to  `f/f'`. However, this uses `f' ~ fp = (fx - f(x-fx))/fx` (a Steffensen step). In
this implementation, `Order1B`, when `fx` is too big, a single secant step of `f`
is used.

The *asymptotic* error, `eᵢ = xᵢ - α`, is given by
`eᵢ₊₂ = 1/2⋅G''/G'⋅ eᵢ⋅eᵢ₊₁ + (1/6⋅G'''/G' - (1/2⋅G''/G'))^2⋅eᵢ⋅eᵢ₊₁⋅(eᵢ+eᵢ₊₁)`.

"""
struct King <: AbstractSecant end
struct Order1B <: AbstractSecant end

struct KingState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    G0::S
end


function init_state(M::Union{King, Order1B}, F, x₀, x₁::T, fx₀, fx₁::S) where {T,S}
    fₛ₀ = F(x₀ - fx₀*oneunit(T)/oneunit(S))
    G₀ = -fx₀^2 / (fₛ₀ - fx₀)
    KingState(x₁, x₀, fx₁, fx₀, G₀)
end
initial_fncalls(::Union{King, Order1B}) = 3

function update_state(method::Order1B, F, o::KingState{T,S}, options, l=NullTracks())  where {T, S}
    if do_guarded_step(Order1B(), o)
        state, flag = update_state(Order1(), F, o, options, l)

        x0,fx0 = state.xn0, state.fxn0 # clunky! Need to update G₀ after Order1() step
        fₛ = F(x0 - fx0*oneunit(T)/oneunit(S))
        incfn(l)
        G₀ = -fx0^2 / (fₛ - fx0)
        @set! state.G0 = G₀

        return (state, flag)
    else
        update_state(King(), F, o, options, l)
    end
end


function update_state(method::King, F,
                      o::KingState{T,S}, options, l=NullTracks())  where {T, S}


    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1
    G₀ = o.G0
    fₛ₁ = F(x1 - fx1*oneunit(T)/oneunit(S))
    incfn(l)
    G₁ = -fx1^2 / (fₛ₁ - fx1)

    m = (x1-x0) / (G₁-G₀) # approximate value of `m`, the multiplicity

    if abs(m) <= 1e-2 * oneunit(m)
        log_message(l, "Estimate for multiplicity has issues. ")
        return (o, true)
    end

    Δ = G₁ * (x1 - x0) / (G₁ - G₀)

    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1 = x1, x1 - Δ
    fx0, fx1 = fx1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1
    @set! o.G0 = G₁


    return o, false

end

##################################################

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


Example
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
struct Esser <: AbstractSecant end
struct Order2B <: AbstractSecant end

function update_state(method::Order2B, fs, o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}
     update_state_guarded(method, Secant(), Esser(), fs, o, options, l)
end

function update_state(method::Esser, F,
                      o::AbstractUnivariateZeroState{T,S}, options, l=NullTracks())  where {T, S}


    x1, fx1 = o.xn1, o.fxn1

    f0 = fx1

    f1  = F(x1 + f0 * oneunit(T) / oneunit(S))
    f_1 = F(x1 - f0 * oneunit(T) / oneunit(S))
    incfn(l, 2)

    # h = f0
    # r1 = f/f' ~ f/f[x+h,x-h]
    # r2 = f'/f'' ~ f[x+h, x-h]/f[x-h,x,x+h]
    r1 = f0 * 2*f0 / (f1 - f_1) * oneunit(T) / oneunit(S)
    r2 = (f1 - f_1)/(f1 - 2*f0 + f_1) * f0/2 * oneunit(T) / oneunit(S)

    k = r2/(r2-r1)  # ~ m

    if abs(k) <= 1e-2 * oneunit(k)
        log_message(l, "Estimate for multiplicity had issues. ")
        return o, true
    end

    delta = k * r1

    if isissue(delta)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1 = x1, x1 - delta
    fx0, fx1 = fx1, F(x1)
    incfn(l)



    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false

end


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
struct Schroder <: AbstractHalleyLikeMethod
end
const Schroeder = Schroder # either spelling
const Schröder = Schroder

## Shroder-like methods
function init_state(M::Schroder, F::Callable_Function, x)
    x₀ = float(first(x))
    fx₀, (Δ, ΔΔ) = F(x₀)
    x₁ = x₀ - ΔΔ/(ΔΔ - Δ) * Δ # m*r1
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end



function update_state(method::Schroder, F, o::AbstractUnivariateZeroState{T,S},
                      options::UnivariateZeroOptions, l=NullTracks()) where {T,S}
    xn = o.xn1
    fxn = o.fxn1
    r1, r2 = o.Δ, o.ΔΔ

    Δ =  r2 / (r2 - r1) * r1  # m * r1

    if isissue(Δ)
        log_message(l, "Issue with increment")
        return o, true
    end

    xn1::T = xn - Δ

    fxn1::S, (r1::T, r2::T) = F(xn1)
    incfn(l,3)

    @set! o.xn0 = xn
    @set! o.xn1 = xn1
    @set! o.fxn0 = fxn
    @set! o.fxn1 = fxn1
    @set! o.Δ = r1
    @set! o.ΔΔ = r2

    return o, false
end


## Thukral 3,4,5 (2) is Schroder
"""
    AbstractThukralBMethod

Abstract type for `ThukralXB` methods for `X` being `2`,`3`,`4`, or `5`.

These are a family of methods which are
* efficient (order `X`) for non-simple roots (e.g. `Thukral2B` is the `Schroder` method)
* taking `X+1` function calls per step
* require `X` derivatives. These can be passed as a tuple of functions, `(f, f', f'', …)`, *or* as
a function returning the ratios: `x -> (f(x), f(x)/f'(x), f'(x)/f''(x), …)`.

## Example

```julia
using ForwardDiff
Base.adjoint(f::Function)  = x  -> ForwardDiff.derivative(f, float(x))
f(x) = (exp(x) + x - 2)^6
x0 = 1/4
find_zero((f, f', f''), x0, Roots.Halley())               # 14 iterations; 45 function evaluations
find_zero((f, f', f''), big(x0), Roots.Thukral2B())       #  4 iterations; 15 function evaluations
find_zero((f, f', f'', f'''), big(x0), Roots.Thukral3B()) #  3 iterations; 16 function evaluations
```


## Refrence

*Introduction to a family of Thukral ``k``-order method for finding multiple zeros of nonlinear equations*,
R. Thukral, JOURNAL OF ADVANCES IN MATHEMATICS 13(3):7230-7237, DOI: [10.24297/jam.v13i3.6146](https:doi.org/10.24297/jam.v13i3.6146).
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
    Δs::NTuple{N, T}
    fxn1::S
    fxn0::S
end

function init_state(M::AbstractThukralBMethod, F::Callable_Function, x)
    x₁ = float(first(x))
    fx₁, Δs = F(x₁)
    state = init_state(M, F, nan(x₁), x₁, nan(fx₁), fx₁; Δs=Δs)
end

function init_state(M::AbstractThukralBMethod, F, x₀, x₁::T, fx₀, fx₁; Δs = nothing) where {T}
    ThukralBState(x₁, x₀, NTuple{fn_argout(M)-1,T}(Δs), fx₁, fx₀)
end

function update_state(M::AbstractThukralBMethod, F, o::AbstractUnivariateZeroState{T,S},
                      options::UnivariateZeroOptions, l=NullTracks()) where {T,S}

    x₀ = o.xn1

    Δ = compute_thukral_Δ(M, o)

    isissue(Δ) && return (o, true)
    x₁ = x₀ - Δ
    fx₁, Δs = F(x₁)
    incfn(l, fn_argout(M))

    @set! o.xn0 = x₀
    @set! o.fxn0 = o.fxn1
    @set! o.Δs =  NTuple{fn_argout(M)-1,T}(Δs)
    @set! o.xn1 = x₁
    @set! o.fxn1 = fx₁

    return (o, false)
end


function compute_thukral_Δ(M::Thukral2B, o)
    r₁, r₂ = o.Δs
    t₁, t₂ = 1/r₁, 1/r₂
    Δ = one(o.xn1)
    Δ /= (t₁ - t₂)
    Δ
end


function compute_thukral_Δ(M::Thukral3B, o)
    r₁, r₂, r₃ = o.Δs
    t₁, t₂, t₃ = 1/r₁, 1/r₂, 1/r₃
    Δ = (2t₁ - 2t₂)
    Δ /= (2t₁^2 - 3t₁*t₂ + t₂*t₃)
    Δ
end

function compute_thukral_Δ(M::Thukral4B, o)
    r₁, r₂, r₃, r₄ = o.Δs
    t₁, t₂, t₃, t₄ = 1/r₁, 1/r₂, 1/r₃, 1/r₄
    Δ = 6t₁^2 - 9t₁*t₂ +3t₂*t₃
    Δ /= 6t₁^3 - 12 * t₁^2*t₂ + 4t₁*t₂*t₃ - t₂*t₃*t₄ + 3*t₁*t₂^2
    Δ
end

function compute_thukral_Δ(M::Thukral5B, o)
    r₁, r₂, r₃, r₄, r₅ = o.Δs
    t₁, t₂, t₃, t₄, t₅ = 1/r₁, 1/r₂, 1/r₃, 1/r₄, 1/r₅
    Δ = 24*t₁^3 - 48t₁^2*t₂  + 16*t₁*t₂*t₃ - 4*t₂*t₃*t₄ + 12t₁*t₂^2
    Δ /= 24t₁^4  - 60t₁^3*t₂ + 20*t₁^2*t₂*t₃ - 5*t₁*t₂*t₃*t₄ + 30t₁^2*t₂^2 - 10*t₁*t₂^2*t₃ + t₂*t₃*t₄*t₅
    Δ
end
