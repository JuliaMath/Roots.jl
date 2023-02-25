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
struct Secant <: AbstractSecantMethod end
const Order1 = Secant
const Orderφ = Secant

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

    x0::T, x1::T = xn1, xn1 - Δ
    fx0::S, fx1::S = fxn1, F(x1)
    incfn(l)

    o = _set(o, (x1, fx1), (x0, fx0))

    return o, false
end
