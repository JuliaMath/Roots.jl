struct FalsePosition{R} <: AbstractSecantMethod end

"""

    FalsePosition([galadino_factor])

Use the [false
position](https://en.wikipedia.org/wiki/False_position_method) method
to find a zero for the function `f` within the bracketing interval
`[a,b]`.

The false position method is a modified bisection method, where the
midpoint between `[aₖ, bₖ]` is chosen to be the intersection point
of the secant line with the ``x`` axis, and not the average between the
two values.

To speed up convergence for concave functions, this algorithm
implements the ``12`` reduction factors of Galdino (*A family of regula
falsi root-finding methods*). These are specified by number, as in
`FalsePosition(2)` or by one of three names `FalsePosition(:pegasus)`,
`FalsePosition(:illinois)`, or `FalsePosition(:anderson_bjork)` (the
default). The default choice has generally better performance than the
others, though there are exceptions.

For some problems, the number of function calls can be greater than
for the `Bisection` method, but generally this algorithm will make
fewer function calls.

Examples
```
find_zero(x -> x^5 - x - 1, (-2, 2), FalsePosition())
```
"""
FalsePosition
FalsePosition(x=:anderson_bjork) = FalsePosition{x}()

# 12 is tough; needs more evaluations
function default_tolerances(::FalsePosition{12}, ::Type{T}, ::Type{S}) where {T,S}
    xatol = eps(real(T)) * oneunit(real(T))
    xrtol = eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4 * eps(real(float(S))) * one(real(S))
    maxiters = 250
    strict = false
    (xatol, xrtol, atol, rtol, maxiters, strict)
end

init_state(M::FalsePosition, F, x₀, x₁, fx₀, fx₁) =
    init_state(Bisection(), F, x₀, x₁, fx₀, fx₁)

function update_state(
    method::FalsePosition,
    fs,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    a, b = o.xn0, o.xn1
    fa, fb = o.fxn0, o.fxn1
    lambda = fb / (fb - fa)

    ϵ = √eps(T) / 100 # some engineering to avoid short moves; still fails on some
    ϵ ≤ lambda ≤ 1 - ϵ || (lambda = 1 / 2)
    x::T = b - lambda * (b - a)
    fx::S = fs(x)
    incfn(l)

    iszero(fx) && return (_set(o, (x, fx)), true)

    if sign(fx) * sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino_reduction(method, fa, fb, fx)
    end
    b, fb = x, fx

    o = _set(o, (b, fb), (a, fa))

    return (o, false)
end

# the 12 reduction factors offered by Galdino
# In RootsTesting.jl, we can see :12 has many more failures.
galdino = Dict{Union{Int,Symbol},Function}(
    :1 => (fa, fb, fx) -> fa * fb / (fb + fx),
    :2 => (fa, fb, fx) -> (fa - fb) / 2,
    :3 => (fa, fb, fx) -> (fa - fx) / (2 + fx / fb),
    :4 => (fa, fb, fx) -> (fa - fx) / (1 + fx / fb)^2,
    :5 => (fa, fb, fx) -> (fa - fx) / (1.5 + fx / fb)^2,
    :6 => (fa, fb, fx) -> (fa - fx) / (2 + fx / fb)^2,
    :7 => (fa, fb, fx) -> (fa + fx) / (2 + fx / fb)^2,
    :8 => (fa, fb, fx) -> fa / 2,
    :9 => (fa, fb, fx) -> fa / (1 + fx / fb)^2,
    :10 => (fa, fb, fx) -> (fa - fx) / 4,
    :11 => (fa, fb, fx) -> fx * fa / (fb + fx),
    :12 => (fa, fb, fx) -> (fa * (1 - fx / fb > 0 ? 1 - fx / fb : 1 / 2)),
)

# give common names
for (nm, i) in [(:pegasus, 1), (:illinois, 8), (:anderson_bjork, 12)]
    galdino[nm] = galdino[i]
end

# from Chris Elrod; https://raw.githubusercontent.com/chriselrod/AsymptoticPosteriors.jl/master/src/false_position.jl
@generated function galdino_reduction(methods::FalsePosition{R}, fa, fb, fx) where {R}
    f = galdino[R]
    quote
        $Expr(:meta, :inline)
        $f(fa, fb, fx)
    end
end
