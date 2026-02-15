# Thanks to Proektsoftbg for opening https://github.com/JuliaMath/Roots.jl/issues/487
# and for providing the basic code!!!

"""
    Roots.modAB()

Implementation of Modified Anderson-Bjork’s method for solving non-linear equations in structural mechanics by N Ganchovski and A Traykov.

This method is a "hybrid" method which chooses between Anderson-Bjork or Bisection.

It is intended to fast and stable enough for equations in structural mechanics.

Andeson-Bjork is a false-position method with adjustments to avoid endpoints that don't move. It has super-linear convergence with a factor of 1.7. This method starts with (simple) bisection and switches to Anderson-Bjork when the difference between the value at the midpoint and the value of the midline are smalle enough.

## Examples

```
```

## Reference

N Ganchovski and A Traykov 2023 IOP Conf. Ser.: Mater. Sci. Eng. 1276 012010

DOI 10.1088/1757-899X/1276/1/012010

[https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010](https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010)

## Hat tip:
Thanks much to `Proektsoftbg` for suggesting the method and providing an implementation (https://github.com/Proektsoftbg/Numerical/tree/main/Numerical-Julia)

"""
struct ModAB <: AbstractBracketingMethod end

# init state
struct ModABState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
    cnt::Int
    bisection::Bool
    side::Symbol
end

function init_state(M::ModAB, F::Callable_Function, x)
    x₀, x₁ = extrema(x₀x₁(x))
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))

    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

# compute fx₁, Δ
function init_state(::ModAB, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    cnt = 0
    bisection = true
    side = :nothing
    assert_bracket(fx₀, fx₁)
    ModABState(promote(x₁, x₀)..., promote(fx₁, fx₀)..., cnt, bisection, side)
end



initial_fncalls(M::ModAB) = 2

function update_state(
    ::ModAB,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}

    x1, x2 = o.xn0, o.xn1
    y1, y2 = o.fxn0, o.fxn1

    # The results show that there is no significant difference (about one iteration)
    # for k (0.1, 0.9). The optimal value is located between 0.1 – 0.4.
    κ = one(x1) / 4

    N = -log2(eps(T))/2 + 1

    cnt::Int = o.cnt + 1
    bisection::Bool = o.bisection
    side::Symbol = o.side

    if bisection
        # take bisection step
        x3 = x1/2 + x2/2 #__middle(x1, x2)
        y3::S = F(x3)
        incfn(l)
        ym = y1/2 + y2/2
        if abs(ym - y3) < κ * (abs(ym) + abs(y3))
            bisection = false
        end
    else
        # take false position step
        x3 = (x1*y2 - y1*x2) / (y2 - y1)
        y3 = F(x3)
        incfn(l)
    end


    δ, ϵ = options.xabstol, options.xreltol
    if iszero(y3) ||
        (side == :left && isapprox(x3, x1; atol=δ, rtol=ϵ)) ||
        (side == :right && isapprox(x3, x2; atol=δ, rtol=ϵ))
        @reset o.xn0 = x2
        @reset o.xn1 = x3
        @reset o.fxn0 = y2
        @reset o.fxn1 = y3

        return (o, true)
    end


    if side == :left
        m = 1 - y3/y1
        y2 = m ≤ 0 ? y2/2 : y2*m)
    elseif side == :right
        m = 1 - y3/y2
        y1 = m ≤ 0 ? y1/2 : y1*m)
    end

    if sign(y1) == sign(y3)
        !bisection && (side = :left)
        x1, y1 = x3, y3
    else
        !bisection && (side = :right)
        x2, y2 = x3, y3
    end

    if rem(cnt, N) == 0
        bisection == true # restart
        side = :nothing
    end


    @reset o.xn0 = x1
    @reset o.xn1 = x2
    @reset o.fxn0 = y1
    @reset o.fxn1 = y2
    @reset o.cnt = cnt
    @reset o.bisection = bisection
    @reset o.side = side

    return o, false

end




#=

    x1, x2 = min(left, right), max(left, right)
    y1 = f(x1) - target
    abs(y1) <= precision && return x1

    y2 = f(x2) - target
    abs(y2) <= precision && return x2

    n_max = -Int(floor(log2(precision) / 2.0)) + 1
    eps1 = precision / 4
    eps = precision * (x2 - x1) / 2.0
    if abs(target) > 1
        eps1 *= target
    end

    side = 0
    ans = x1
    bisection = true
    k = 0.25

    for i in 1:200
        if bisection
            x3 = (x1 + x2) / 2.0
            y3 = f(x3) - target
            ym = (y1 + y2) / 2.0
            if abs(ym - y3) < k * (abs(y3) + abs(ym))
                bisection = false
            end
        else
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
            if x3 < x1 - eps || x3 > x2 + eps
                return NaN
            end
            y3 = f(x3) - target
        end

        if abs(y3) < eps1 || abs(x3 - ans) < eps
            if x1 > x2
                return side == 1 ? x2 : x1
            end
            return clamp(x3, x1, x2)
        end

        ans = x3
        if sign(y1) == sign(y3)
            if side == 1
                m = 1 - y3 / y1
                if m <= 0
                    y2 /= 2
                else
                    y2 *= m
                end
            elseif !bisection
                side = 1
            end
            x1 = x3
            y1 = y3
        else
            if side == -1
                m = 1 - y3 / y2
                if m <= 0
                    y1 /= 2
                else
                    y1 *= m
                end
            elseif !bisection
                side = -1
            end
            x2 = x3
            y2 = y3
        end

        if i % n_max == 0
            bisection = true
        end
    end
    return ans
end
=#
