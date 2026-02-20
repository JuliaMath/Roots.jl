# Thanks to Proektsoftbg for opening https://github.com/JuliaMath/Roots.jl/issues/487
# and for providing the basic code!!!

"""
    Roots.modAB()

Implementation of a modified Anderson-Bjork method by N Ganchovski and A Traykov.

The method is advertised as intended to be fast and stable enough for solving non-linear equations in structural mechanics.

This method is a "hybrid" method which chooses between Anderson-Bjork or simple bisection, starting with bisection and moving to Anderson-Bjork when a condition is met. The condition used is a measure of the relative deviation of the function from a straight line at the middle point between the bracketing values. If too many steps are taken, the  more robust bisection method is used to ensure convergence in a bounded number of steps.

Andeson-Bjork is a false-position method with adjustments to avoid endpoints that don't move. It has super-linear convergence with a factor of `≈ 1.7`.

## Examples

```
julia> mutable struct Cnt
           cnt::Int
           f
           Cnt(f) = new(0, f)
       end

julia> (f::Cnt)(x) = (f.cnt += 1; f.f(x))

julia> f(x) = (sin(x) - x / 4)^3; (a,b) = (2,4); # multiplicity at root

julia> F = Cnt(f); x = find_zero(F, (a,b), Roots.Bisection()); (x, F.cnt)
(2.4745767873698292, 54)

julia> F = Cnt(f); x = find_zero(F, (a,b), Roots.Ridders()); (x, F.cnt)
(2.474576787369829, 85)

julia> F = Cnt(f); x = find_zero(F, (a,b), Roots.FalsePosition(:anderson_bjork)); (x, F.cnt)
(2.4745666938402757, 40)

julia> F = Cnt(f); x = find_zero(F, (a,b), Roots.ModAB()); (x, F.cnt)
(2.474578857421875, 18)
```

## Notes

* Unlike some other bracketing methods, this implementation stops on a small residual (small `f(xₙ)`, by default with only an absolute tolerance) in addition to a small bracket (small `|xₙ - xₙ₋₁|`). All of `xabstol`, `xreltol`, `abstol`, and `reltol` may be specified. The method, unlike `Bisection`, does not expect that `nextfloat(a) ≥ b` for the last identified bracketing interval `[a,b]`.

* Over `Float64`, this method should handle infinite intervals (likely falling back to the method of `Bisection`).

* This method seems to perform quite efficiently when the root has multiplicities, as in the example.

* That the value found in the example by `ModAB` and `Bisection` agree only through the first 5 decimal points is due to the `ModAB` algorithm stopping on small `f(xₙ)` values, as `Bisection` iterates up to the last floating point bit unless it finds an exact numeric zero.

## Reference

N Ganchovski and A Traykov 2023 IOP Conf. Ser.: Mater. Sci. Eng. 1276 012010

DOI 10.1088/1757-899X/1276/1/012010

[https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010](https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010)


A new high order method of regula falsi type for computing a root of an equation; Ned Anderson & Åke Björck. [https://link.springer.com/article/10.1007/BF01951936](https://link.springer.com/article/10.1007/BF01951936).

## Hat tip:

Thanks much to `Proektsoftbg` for suggesting the method and providing an implementation (https://github.com/Proektsoftbg/Numerical/tree/main/Numerical-Julia)

"""
struct ModAB <: AbstractSecantMethod end

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
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))

    state = init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function init_state(::ModAB, F, x₀::T, x₁::T, fx₀, fx₁) where {T}
    cnt = 0
    bisection = true
    side = :nothing
    assert_bracket(fx₀, fx₁)
    ModABState(promote(x₁, x₀)..., promote(fx₁, fx₀)..., cnt, bisection, side)
end

initial_fncalls(M::ModAB) = 2

# adjust to stop on smallish Δx
# xatol is small, but allows convergence of extreme case:
# find_zero(x ->  x < eps(0.0) ? -1.0 : 1.0, (-Inf, Inf), Roots.ModAB())
# maxiters for T ∈ (Float16, Float32, Float64) need only be N + X
# where X is 16,32,or 64, N is -Int(log2(eps(T))) ÷ 2 + 1
function default_tolerances(::ModAB,::Type{T},::Type{S}) where {T,S}
    xatol = 2 * eps(zero(T)) * oneunit(real(T)) # not quite 0
    xrtol = eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = zero(real(S))
    maxiters = 3 * (-Int(log2(eps(T))) + 1)
    strict = true
    (xatol, xrtol, atol, rtol, maxiters, strict)
end


function update_state(
    ::ModAB,
    F,
    o::AbstractUnivariateZeroState{T,S},
    options,
    l=NullTracks(),
) where {T,S}

    x1, x2 = o.xn0, o.xn1
    y1, y2 = o.fxn0, o.fxn1

    if x1 > x2 # we store x3 as x2; fix order
        x1, x2 = x2, x1
        y1, y2 = y2, y1
    end


    # The results show that there is no significant difference (about one iteration)
    # for k (0.1, 0.9). The optimal value is located between 0.1 – 0.4.
    κ = one(x1) / 4


    cnt::Int = o.cnt + 1
    bisection::Bool = o.bisection
    side::Symbol = o.side
    gaveup = false
    N = -Int(log2(eps(T))) ÷ 2 + 1
    cnt == N && @show :giveup, N
    # find x3, y3

    if cnt > N
        gaveup = true # finish with Bisection() method
        if sign(x1) * sign(x2) < 0
            x3 = zero(x1)
        else
            x3 = __middle(x1, x2)
        end
        y3 = first(F(x3))
        incfn(l)
    elseif bisection
        # take bisection step
        x3 = x1/2 + x2/2
        y3 = first(F(x3))
        incfn(l)

        # continue with simple bisection?
        ym = y1/2 + y2/2
        if abs(ym - y3) < κ * (abs(ym) + abs(y3))
            bisection = false
        end
    else
        # take false position step
        x3 = (x1*y2 - y1*x2) / (y2 - y1)
        y3 = first(F(x3))
        incfn(l)
    end

    # anderson-bjork adjustment to y
    if sign(y1) == sign(y3)
        if side == :right && !gaveup
            m = 1 - y3 / y1
            y2 = m ≤ 0 ? y2/2 : y2 * m
        elseif !bisection
            side = :right
        end
        # use x3,x2
    else
        if side == :left && !gaveup
            m = 1 - y3 / y2
            y1 = m ≤ 0 ? y1/2 : y1*m
        elseif !bisection
            side = :left
        end
        #use x1 x2
        x2, y2 = x1, y1
    end

    @reset o.xn0 = x2
    @reset o.xn1 = x3
    @reset o.fxn0 = y2
    @reset o.fxn1 = y3
    @reset o.cnt = cnt
    @reset o.bisection = bisection
    @reset o.side = side

    return o, false

end
