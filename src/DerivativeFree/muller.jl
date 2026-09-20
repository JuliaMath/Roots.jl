"""
    Muller()

*Muller’s method* generalizes the secant method, but uses quadratic
interpolation among three points instead of linear interpolation between two.
Solving for the zeros of the quadratic allows the method to find complex
pairs of roots.

Given three previous guesses for the root `xᵢ₋₂`, `xᵢ₋₁`, `xᵢ`, and the values
of the polynomial `f` at those points, the next approximation `xᵢ₊₁` is produced.

Excerpt and the algorithm taken from W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery *Numerical Recipes in C*, Cambridge University Press (2002), p. 371.

Muller's method can return complex zeros if the initial guess is complex.

"""
struct Muller <: AbstractSecantMethod end
initial_fncalls(::Muller) = 3

struct MullerState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn2::T
    xn1::T
    xn0::T
    fxn2::S
    fxn1::S
    fxn0::S
end

function _muller_bootstrap(f::Callable_Function, x::Number)
    a, b = x₀x₁(float(x))
    fa, fb = f.((a, b))
    c = b - fb * (b-a)/(fb-fa)
    fc = f(c)
    (a, b, c, fa, fb, fc)
end
function _muller_bootstrap(f::Callable_Function, x1::Number, x2::Number)
    a, b = x₀x₁((float(x1), float(x2)))
    fa, fb = f.((a, b))
    c = b - fb * (b-a)/(fb-fa)
    fc = f(c)
    (a, b, c, fa, fb, fc)
end
function _muller_bootstrap(f::Callable_Function, x1::Number, x2::Number, x3)
    a, b, c = float.((x1, x2, x3))
    fa, fb, fc = f(a), f(b), f(c)
    (a, b, c, fa, fb, fc)
end

function init_state(M::Muller, F::Callable_Function, x)
    # can specify 1, 2 or 3 initial points
    x0, x1, x2, fx0, fx1, fx2 = _muller_bootstrap(F, x...)
    state = MullerState(x2, x1, x0, fx2, fx1, fx0)
end

function update_state(
    M::Muller,
    F,
    o::MullerState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    a, b, c = o.xn0, o.xn1, o.xn2
    fa, fb, fc = o.fxn0, o.fxn1, o.fxn2

    q = (c - b) / (b - a)
    q² = q^2
    q1 = q + one(q)

    A = q * fc - q * q1 * fb + q² * fa
    B = (q1 + q) * fc - q1^2 * fb + q² * fa
    C = q1 * fc

    Δ = B^2 - 4A * C
    if typeof(Δ) <: Real &&
       Δ < 0 &&
       throw(
           DomainError(
               Δ,
               "Discriminant is negative and the function most likely has complex roots. You might use the `Roots.Muller()` method with complex input.",
           ),
       )
    end

    Δ = √Δ
    d⁺ = B + Δ
    d⁻ = B - Δ
    den = abs(d⁺) > abs(d⁻) ? d⁺ : d⁻
    inc = (c - b) * 2C / den

    if isissue(inc)
        log_message(l, "Increment `Δx` has issues. ")
        return o, true
    end

    x = c - inc
    fx = first(F(x))

    # check for oddity?

    incfn(l)
    @reset o.xn0 = b
    @reset o.xn1 = c
    @reset o.xn2 = x
    @reset o.fxn0 = fb
    @reset o.fxn1 = fc
    @reset o.fxn2 = fx

    return (o, false)
end
