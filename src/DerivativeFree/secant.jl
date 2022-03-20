"""
    Secant()
    Order1()
    Orderφ()


The `Order1()` method is an alias for `Secant`. It specifies the
[secant method](https://en.wikipedia.org/wiki/Secant_method).
This method keeps two values in its state, `xₙ` and `xₙ₋₁`. The
updated point is the intersection point of ``x`` axis with the secant line
formed from the two points. The secant method uses ``1`` function
evaluation per step and has order `φ≈ (1+sqrt(5))/2`.

The error, `eᵢ = xᵢ - α`, satisfies
`eᵢ₊₂ = f[xᵢ₊₁,xᵢ,α] / f[xᵢ₊₁,xᵢ] * (xᵢ₊₁-α) * (xᵢ - α)`.

"""
struct Secant <: AbstractSecant end
const Order1 = Secant
const Orderφ = Secant

initial_fncalls(::AbstractSecant) = 2

# init_state(M,F,x) --> call init_state(M,F,x₀,x₁,fx₀, fx₁)
function init_state(M::AbstractSecant, F::Callable_Function, x)
    x₀, x₁ = x₀x₁(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# initialize from xs, fxs
function init_state(::AbstractSecant, F, x₀, x₁, fx₀, fx₁)
    UnivariateZeroState(x₁, x₀, fx₁, fx₀)
end

function update_state(
    ::Order1,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    xn0, xn1 = o.xn0, o.xn1
    fxn0, fxn1 = o.fxn0, o.fxn1
    Δ = fxn1 * (xn1 - xn0) / (fxn1 - fxn0)

    if isissue(Δ)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1::T = xn1, xn1 - Δ
    fx0, fx1 = fxn1, F(x1)
    incfn(l)

    @set! o.xn0 = x0
    @set! o.xn1 = x1
    @set! o.fxn0 = fx0
    @set! o.fxn1 = fx1

    return o, false
end
