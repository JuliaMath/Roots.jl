# 1 step taken in set up
function log_step(l::Tracks, M::AbstractDerivativeMethod, state; init=false)
    init && push!(l.xfₛ, (state.xn0, state.fxn0))
    push!(l.xfₛ, (state.xn1, state.fxn1))
    init && log_iteration(l, 1)
    !init && log_iteration(l, 1)
    nothing
end

"""

    Roots.Newton()

Implements Newton's [method](https://en.wikipedia.org/wiki/Newton%27s_method):
`xᵢ₊₁ =  xᵢ - f(xᵢ)/f'(xᵢ)`.  This is a quadratically convergent method requiring
one derivative and two function calls per step.

## Examples

```jldoctest with_derivative
julia> using Roots

julia> find_zero((sin,cos), 3.0, Roots.Newton()) ≈ π
true
```

If function evaluations are expensive one can pass in a function which returns (f, f/f') as follows

```jldoctest with_derivative
julia> find_zero(x -> (sin(x), sin(x)/cos(x)), 3.0, Roots.Newton()) ≈ π
true
```

This can be advantageous if the derivative is easily computed from the
value of f, but otherwise would be expensive to compute.

----

The error, `eᵢ = xᵢ - α`, can be expressed as `eᵢ₊₁ =
f[xᵢ,xᵢ,α]/(2f[xᵢ,xᵢ])eᵢ²` (Sidi, Unified treatment of regula falsi,
Newton-Raphson, secant, and Steffensen methods for nonlinear
equations).

"""
struct Newton <: AbstractNewtonLikeMethod end

fn_argout(::AbstractNewtonLikeMethod) = 2

# we store x0,x1,fx0,fx1 **and** Δ = fx1/f'(x1)
struct NewtonState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    Δ::T
    fxn1::S
    fxn0::S
end

function init_state(M::Newton, F::Callable_Function, x)
    x₀ = float(first(x))
    T = eltype(x₀)
    fx₀, Δ::T = F(x₀)
    x₁::T = x₀ - Δ
    state = init_state(M, F, x₀, x₁, fx₀, fx₀)
end

# compute fx₁, Δ
function init_state(::Newton, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    fx₁, Δ::T = F(x₁)
    NewtonState(x₁, x₀, Δ, fx₁, fx₀)
end

initial_fncalls(M::Newton) = 2

function update_state(
    M::Newton,
    F,
    o::NewtonState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ::T = o.Δ

    if isissue(Δ)
        log_message(l, "Issue with `f/f′'")
        return o, true
    end

    xn0, xn1::T = xn1, xn1 - Δ
    fxn0 = fxn1
    fxn1::S, Δ = F(xn1)
    incfn(l, 2)

    @set! o.xn0 = xn0
    @set! o.xn1 = xn1
    @set! o.Δ = Δ
    @set! o.fxn0 = fxn0
    @set! o.fxn1 = fxn1

    return o, false
end
