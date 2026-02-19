# Thanks to Proektsoftbg for opening https://github.com/JuliaMath/Roots.jl/issues/487
# and for providing the basic code!!!

"""
    Roots.modAB()

Implementation of modified Anderson-Bjork method for solving non-linear equations in structural mechanics by N Ganchovski and A Traykov.

This method is a "hybrid" method which chooses between Anderson-Bjork or Bisection.

It is intended to fast and stable enough for equations in structural mechanics.

Andeson-Bjork is a false-position method with adjustments to avoid endpoints that don't move. It has super-linear convergence with a factor of 1.7. This method starts with (simple) bisection and switches to Anderson-Bjork when the difference between the value at the midpoint and the value of the midline are small enough.

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

## Reference

N Ganchovski and A Traykov 2023 IOP Conf. Ser.: Mater. Sci. Eng. 1276 012010

DOI 10.1088/1757-899X/1276/1/012010

[https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010](https://iopscience.iop.org/article/10.1088/1757-899X/1276/1/012010)

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
function default_tolerances(::ModAB,::Type{T},::Type{S}) where {T,S}
    xatol = 4 * eps(real(T)) * oneunit(real(T))
    xrtol = 4 * eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = zero(real(S))
    maxiters = 100
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

    # find x3, y3
    if bisection
        # take bisection step
        x3 = x1/2 + x2/2
        y3::S = F(x3)
        incfn(l)

        # continue with bisection?
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

    # anderson-bjork adjustment to y
    if sign(y1) == sign(y3)
        if side == :right
            m = 1 - y3 / y1
            y2 = m ≤ 0 ? y2/2 : y2 * m
        elseif !bisection
            side = :right
        end
        # use x3,x2
    else
        if side == :left
            m = 1 - y3 / y2
            y1 = m ≤ 0 ? y1/2 : y1*m
        elseif !bisection
            side = :left
        end
        #use x1 x2
        x2, y2 = x1, y1
    end

    # restart bisection?
    N = -Int(log2(eps(T))/2) + 1
    if rem(cnt, N) == 0
        bisection == true
        side = :nothing
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
