"""
    Steffensen()
    Order2()


The quadratically converging
[Steffensen](https://en.wikipedia.org/wiki/Steffensen's_method#Simple_description)
method is used for the derivative-free `Order2()` algorithm. Unlike
the quadratically converging Newton's method, no derivative is
necessary, though like Newton's method, two function calls per step
are. Steffensen's algorithm is more sensitive than Newton's method to
poor initial guesses when `f(x)` is large, due to how `f'(x)` is
approximated. The `Order2` method replaces a Steffensen step with a secant
step when `f(x)` is large.

The error, `eᵢ - α`, satisfies
`eᵢ₊₁ = f[xᵢ, xᵢ+fᵢ, α] / f[xᵢ,xᵢ+fᵢ] ⋅ (1 - f[xᵢ,α] ⋅ eᵢ²`
"""
struct Steffensen <: AbstractSecantMethod end
struct Order2 <: AbstractSecantMethod end

function update_state(
    M::Order2,
    fs,
    o::UnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    update_state_guarded(M, Secant(), Steffensen(), fs, o, options, l)
end

function update_state(
    ::Steffensen,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    x0, x1 = o.xn0, o.xn1
    fx0, fx1 = o.fxn0, o.fxn1

    sgn = sign((fx1 - fx0) / (x1 - x0))
    x2 = x1 - sgn * fx1 / oneunit(S) * oneunit(T)

    f0 = fx1
    f1::S = F(x2)
    incfn(l, 1)

    delta = -sgn * f0 * f0 / (f1 - f0) * oneunit(T) / oneunit(S)

    if isissue(delta)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x0, x1::T = x1, x1 - delta
    fx0, fx1::S = fx1, F(x1)
    incfn(l)

    o = _set(o, (x1, fx1), (x0, fx0))

    return o, false
end
