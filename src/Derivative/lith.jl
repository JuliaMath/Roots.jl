"""
    LithBoonkkampIJzerman{S,D} <: AbstractNewtonLikeMethod
    LithBoonkkampIJzerman(S,D)

A family of different methods that includes the secant method and Newton's method.

Specifies a linear multistep solver with `S` steps and `D` derivatives following [Lith, Boonkkamp, and
IJzerman](https://doi.org/10.1016/j.amc.2017.09.003).

# Extended help

## Examples

```jldoctest lith
julia> using Roots

julia> find_zero(sin, 3, Roots.LithBoonkkampIJzerman(2,0)) ≈ π # the secant method
true

julia> find_zero((sin,cos), 3, Roots.LithBoonkkampIJzerman(1,1)) ≈ π # Newton's method
true

julia> find_zero((sin,cos), 3, Roots.LithBoonkkampIJzerman(3,1)) ≈ π # Faster convergence rate
true

julia> find_zero((sin,cos, x->-sin(x)), 3, Roots.LithBoonkkampIJzerman(1,2)) ≈ π # Halley-like method
true
```

The method can be more robust to the initial condition. This example is from the paper (p13). Newton's method (the `S=1`, `D=1` case) fails if `|x₀| ≥ 1.089` but methods with more memory succeed.

```jldoctest lith
julia> fx =  ZeroProblem((tanh,x->sech(x)^2), 1.239); # zero at 0.0

julia> solve(fx, Roots.LithBoonkkampIJzerman(1,1)) |> isnan# Newton, NaN
true

julia> solve(fx, Roots.LithBoonkkampIJzerman(2,1)) |> abs |> <(eps())
true

julia> solve(fx, Roots.LithBoonkkampIJzerman(3,1)) |> abs |> <(eps())
true
```

Multiple derivatives can be constructed automatically using automatic differentiation. For example,

```jldoctest lith
julia> using ForwardDiff

julia> function δ(f, n::Int=1)
           n <= 0 && return f
           n == 1 && return x -> ForwardDiff.derivative(f,float(x))
           δ(δ(f,1),n-1)
       end;

julia> fs(f,n) = ntuple(i -> δ(f,i-1), Val(n+1));

julia> f(x) = cbrt(x) * exp(-x^2); # cf. Table 6 in paper, α = 0

julia> fx = ZeroProblem(fs(f,1), 0.1147);

julia> opts = (xatol=2eps(), xrtol=0.0, atol=0.0, rtol=0.0); # converge if |xₙ - xₙ₋₁| <= 2ϵ

julia> solve(fx, Roots.LithBoonkkampIJzerman(1, 1); opts...) |> isnan # NaN -- no convergence
true

julia> solve(fx, Roots.LithBoonkkampIJzerman(2, 1); opts...) |> abs |> <(eps()) # converges
true

julia> fx = ZeroProblem(fs(f,2), 0.06);                       # need better starting point

julia> solve(fx, Roots.LithBoonkkampIJzerman(2, 2); opts...) |> abs |> <(eps()) # converges
true
```

For the case `D=1`, a bracketing method based on this approach is implemented in [`LithBoonkkampIJzermanBracket`](@ref)

## Reference

In [Lith, Boonkkamp, and
IJzerman](https://doi.org/10.1016/j.amc.2017.09.003) an analysis is
given of the convergence rates when using linear multistep methods to
solve `0=f(x)` using `f⁻¹(0) = x` when `f` is a sufficiently smooth
linear function. The reformulation, attributed to Grau-Sanchez, finds
a differential equation for `f⁻¹`: `dx/dy = [f⁻¹]′(y) = 1/f′(x) = F` as
`x(0) = x₀ + ∫⁰_y₀ F(x(y)) dy`.

A linear multi-step method is used to solve this equation
numerically.  Let S be the number of memory steps (S= 1,2,...) and D be
the number of derivatives employed, then, with `F(x) = dx/dy`
`x_{n+S} = ∑_{k=0}^{S-1} aₖ x_{n+k} +∑d=1^D ∑_{k=1}^{S-1} aᵈ_{n+k}F⁽ᵈ⁾(x_{n+k})`.
 The `aₖ`s and `aᵈₖ`s are computed each step.

This table is from Tables 1 and 3 of the paper and gives the
convergence rate for simple roots identified therein:

```
s: number of steps remembered
d: number of derivatives uses
s/d  0    1    2    3    4
1    .    2    3    4    5
2    1.62 2.73 3.79 4.82 5.85
3    1.84 2.91 3.95 4.97 5.98
4    1.92 2.97 3.99 4.99 5.996
5    1.97 .    .    .    .
```

That is, more memory leads to a higher convergence rate; more
derivatives leads to a higher convergence rate. However, the
interval about `α`, the zero, where the convergence rate is guaranteed
may get smaller.

!!! note
    For the larger values of `S`, the expressions to compute the next value get quite involved.
    The higher convergence rate is likely only to be of help for finding solutions to high precision.

"""
struct LithBoonkkampIJzerman{S,D} <: AbstractNewtonLikeMethod end
LithBoonkkampIJzerman(s::Int, d::Int) = LithBoonkkampIJzerman{s,d}()
fn_argout(::LithBoonkkampIJzerman{S,D}) where {S,D} = 1 + D

struct LithBoonkkampIJzermanState{S′,D⁺,T,S} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    m::NTuple{S′,T}
    fxn1::S
    fxn0::S
    fm::NTuple{D⁺,NTuple{S′,S}}
end

log_step(l::Tracks, M::LithBoonkkampIJzerman, state; init=false) =
    log_step(l, Secant(), state; init=init)

# return f^(i-1)(x) for i in 0:N-1; not the same as default eval call
function evalf(
    F::Callable_Function{S,T,𝑭,P},
    x,
) where {N,S<:Val{N},T<:Val{true},𝑭,P<:Nothing}
    fi = map(f -> f(x), F.f) #recommended on Slack to not allocate
    R = typeof(float(first(fi)))
    convert(NTuple{N,R}, fi)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {N,S<:Val{N},T<:Val{true},𝑭,P}
    fi = map(f -> f(x, F.p), F.f)
    R = typeof(float(first(fi)))
    convert(NTuple{N,R}, fi)
end

#specializations for N = 1,2,3,4,5,6
## lmm(::Roots.LithBoonkkampIJzerman{1, D}) is defined up unitl D = 6, so specialize those

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{1},T<:Val{false},𝑭,P}
    F(x)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{2},T<:Val{false},𝑭,P}
    f, Δ = F(x)
    f′ = f / Δ
    (f, f′)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{3},T<:Val{false},𝑭,P}
    f, Δ = F(x)
    Δ₁, Δ₂ = Δ
    f′ = f / Δ₁
    f′′ = f′ / Δ₂
    (f, f′, f′′)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{4},T<:Val{false},𝑭,P}
    f, Δ = F(x)
    Δ₁, Δ₂, Δ₃ = Δ
    f′ = f / Δ₁
    f′′ = f′ / Δ₂
    f′′′ = f′′ / Δ₃
    (f, f′, f′′, f′′′)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{5},T<:Val{false},𝑭,P}
    f, Δ = F(x)
    Δ₁, Δ₂, Δ₃, Δ₄ = Δ
    f′ = f / Δ₁
    f′′ = f′ / Δ₂
    f′′′ = f′′ / Δ₃
    f′′′′ = f′′′ / Δ₄
    (f, f′, f′′, f′′′, f′′′′)
end

function evalf(F::Callable_Function{S,T,𝑭,P}, x) where {S<:Val{6},T<:Val{false},𝑭,P}
    f, Δ = F(x)
    Δ₁, Δ₂, Δ₃, Δ₄, Δ₅ = Δ
    f′ = f / Δ₁
    f′′ = f′ / Δ₂
    f′′′ = f′′ / Δ₃
    f′′′′ = f′′′ / Δ₄
    f′′′′′ = f′′′′ / Δ₅
    (f, f′, f′′, f′′′, f′′′′, f′′′′′)
end

#function to obtain just the first value. optimized in case of tuple function
function only_f(F::Callable_Function{S,T,𝑭,P}, x) where {N,S<:Val{N},T<:Val{true},𝑭,P}
    return F.f[1](x, F.p)
end

function only_f(
    F::Callable_Function{S,T,𝑭,P},
    x,
) where {N,S<:Val{N},T<:Val{true},𝑭,P<:Nothing}
    return F.f[1](x)
end

function only_f(F::Callable_Function{S,T,𝑭,P}, x) where {N,S<:Val{N},T<:Val{false},𝑭,P}
    return first(F(x))
end

function init_state(L::LithBoonkkampIJzerman{S,0}, F::Callable_Function, x) where {S}
    x₀, x₁ = x₀x₁(x)
    fx₀, fx₁ = only_f(F, x₀), only_f(F, x₁)
    state = init_state(L, F, x₀, x₁, fx₀, fx₁, nothing)
end

function init_state(L::LithBoonkkampIJzerman{S,D}, F::Callable_Function, x) where {S,D}
    x₀ = float(first(x))
    ys₀ = evalf(F, x₀)
    fx₀ = first(ys₀)
    state = init_state(L, F, nan(x₀), x₀, nan(fx₀), fx₀, ys₀)
end

function init_state(
    L::LithBoonkkampIJzerman{S,D},
    F,
    x₀,
    x₁::R,
    fx₀,
    fx₁::T,
    ys₀,
) where {S,D,R,T}
    xs, ys = init_lith(L, F, x₁, fx₁, x₀, fx₀, ys₀) # [x₀,x₁,…,xₛ₋₁], ...
    # skip unit consideration here, as won't fit within storage of ys
    state = LithBoonkkampIJzermanState{S,D + 1,R,T}(
        xs[end],    # xₙ
        S > 1 ? xs[end - 1] : nan(xs[end]), # xₙ₋₁
        xs,         # all xs
        ys[1][end], # fₙ
        S > 1 ? ys[1][end - 1] : nan(ys[1]), # fₙ₋₁
        ys,        #
    )

    state
end
initial_fncalls(::LithBoonkkampIJzerman{S,D}) where {S,D} = S * (D + 1)

function update_state(
    L::LithBoonkkampIJzerman{S,D},
    F::Callable_Function,
    o::LithBoonkkampIJzermanState{S⁺,D′,R,T},
    options,
    l=NullTracks(),
) where {S,D,S⁺,D′,R,T}
    xs, ys = o.m, o.fm
    xᵢ::R = lmm(L, xs, ys...)
    isissue(o.xn1 - xᵢ) && return (o, true)

    for i in 1:(S - 1)
        @set! xs[i] = xs[i + 1]
    end
    @set! xs[end] = xᵢ

    ysᵢ = evalf(F, xᵢ)
    for i in 0:D
        i′ = i + 1
        for j in 1:(S - 1)
            @set! ys[i′][j] = ys[i′][j + 1]
        end
        yij::T = ysᵢ[i′]
        @set! ys[i′][end] = yij
    end
    incfn(l, 1 + D)

    @set! o.xn0 = o.xn1
    @set! o.xn1 = xᵢ
    @set! o.fxn0 = o.fxn1
    @set! o.fxn1 = ys[1][end]
    @set! o.m = xs
    @set! o.fm = ys

    return (o, false)
end

# manufacture initial xs, ys
# use lower memory terms to boot strap up. Secant uses initial default step
#D=0, generate [x0].x1,...,xs
function init_lith(
    L::LithBoonkkampIJzerman{S,0},
    F::Callable_Function{Si,Tup,𝑭,P},
    x₁::R,
    fx₁::T,
    x₀::R,
    fx₀::T,
    ys₀,
) where {S,Si,Tup,𝑭,P,R,T}
    xs = NTuple{S,R}(ntuple(_ -> one(R), Val(S)))
    yᵢ = NTuple{S,T}(ntuple(_ -> one(T), Val(S)))
    ys = NTuple{1,NTuple{S,T}}((yᵢ,))

    # build up to get S of them
    x0::R = zero(R)
    if isnan(x₀)
        x0 = _default_secant_step(x₁)
        fx0::T = only_f(F, x0)
    else
        x0, fx0 = x₀, fx₀
    end

    @set! xs[1] = x0
    @set! xs[2] = x₁
    @set! ys[1][1] = fx0
    @set! ys[1][2] = fx₁

    # build up xs, ys
    # redundant code, but here to avoid allocations
    S < 3 && return (xs, ys)
    xᵢ = lmm(Val(2), Val(0), xs, ys)
    y1i = only_f(F, xᵢ)
    @set! xs[3] = xᵢ
    @set! ys[1][3] = y1i

    S < 4 && return (xs, ys)
    xᵢ = lmm(Val(3), Val(0), xs, ys)
    y1i = only_f(F, xᵢ)
    @set! xs[4] = xᵢ
    @set! ys[1][4] = y1i

    S < 5 && return (xs, ys)
    xᵢ = lmm(Val(4), Val(0), xs, ys)
    y1i = only_f(F, xᵢ)
    @set! xs[5] = xᵢ
    @set! ys[1][5] = y1i

    S < 6 && return (xs, ys)
    xᵢ = lmm(Val(5), Val(0), xs, ys)
    y1i = only_f(F, xᵢ)
    @set! xs[6] = xᵢ
    @set! ys[1][6] = y1i

    for i in 7:S #3:S
        xᵢ::R = lmm(Val(i - 1), Val(0), xs, ys) # XXX allocates due to runtime i-1
        y1i::T = only_f(F, xᵢ)
        @set! xs[i] = xᵢ
        @set! ys[1][i] = y1i
    end

    xs, ys
end

#D≥1. ignore x₀
function init_lith(
    L::LithBoonkkampIJzerman{S,D},
    F::Callable_Function{Si,Tup,𝑭,P},
    x₁::R,
    fx₁::T,
    x₀::R,
    fx₀::T,
    ys₀,
) where {S,D,Si,Tup,𝑭,P,R,T}
    xs = NTuple{S,R}(ntuple(_ -> one(R), Val(S)))
    yᵢ = NTuple{S,T}(ntuple(_ -> one(T), Val(S)))
    ys = NTuple{D + 1,NTuple{S,T}}(ntuple(_ -> yᵢ, Val(D + 1)))

    @set! xs[1] = x₁
    for j in 1:(D + 1)
        @set! ys[j][1] = ys₀[j]
    end

    # build up to get S of them
    # redundant code, but here to avoid allocations
    S < 2 && return xs, ys
    xᵢ = lmm(Val(1), Val(D), xs, ys)
    @set! xs[2] = xᵢ
    ysᵢ = evalf(F, xᵢ)
    for j in 1:(D + 1)
        @set! ys[j][2] = ysᵢ[j]
    end

    S < 3 && return xs, ys
    xᵢ = lmm(Val(2), Val(D), xs, ys)
    @set! xs[3] = xᵢ
    ysᵢ = evalf(F, xᵢ)
    for j in 1:(D + 1)
        @set! ys[j][3] = ysᵢ[j]
    end

    S < 4 && return xs, ys
    xᵢ = lmm(Val(3), Val(D), xs, ys)
    @set! xs[4] = xᵢ
    ysᵢ = evalf(F, xᵢ)
    for j in 1:(D + 1)
        @set! ys[j][4] = ysᵢ[j]
    end

    for i in 5:S
        xᵢ::R = lmm(Val(i - 1), Val(D), xs, ys) # XXX allocates! clean up
        @set! xs[i] = xᵢ
        ysᵢ = evalf(F, xᵢ)
        for j in 1:(D + 1)
            @set! ys[j][i] = ysᵢ[j]
        end
    end

    return xs, ys
end

"""
    LithBoonkkampIJzermanBracket()

A bracketing method which is a modification of Brent's method due to
[Lith, Boonkkamp, and
IJzerman](https://doi.org/10.1016/j.amc.2017.09.003). The best
possible convergence rate is 2.91.

A function, its derivative, and a bracketing interval need to be specified.

The state includes the 3 points -- a bracket `[a,b]` (`b=xₙ` has
`f(b)` closest to `0`) and `c=xₙ₋₁` -- and the corresponding values
for the function and its derivative at these three points.

The next proposed step is either a `S=2` or `S=3` selection for the
[`LithBoonkkampIJzerman`](@ref) methods with derivative information
included only if it would be of help. The proposed is modified if it
is dithering. The proposed is compared against a bisection step; the
one in the bracket and with the smaller function value is chosen as
the next step.


"""
struct LithBoonkkampIJzermanBracket <: AbstractBracketingMethod end
struct LithBoonkkampIJzermanBracketState{T,S,R} <: AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    c::T
    fxn1::S
    fxn0::S
    fc::S
    fp1::R
    fp0::R
    fpc::R
end

fn_argout(::LithBoonkkampIJzermanBracket) = 2

function init_state(M::LithBoonkkampIJzermanBracket, F::Callable_Function, x)
    x₀, x₁ = adjust_bracket(x)
    fx₀, Δfx₀ = F(x₀)
    fx₁, Δfx₁ = F(x₁)
    a, b, fa, fb, f′a, f′b = x₀, x₁, fx₀, fx₁, fx₀ / Δfx₀, fx₁ / Δfx₁
    if abs(fa) < abs(fb)
        a, b, fa, fb = b, a, fb, fa
    end
    assert_bracket(fa, fb)

    c, fc, f′c = a, fa, f′a

    # skip unit consideration here, as won't fit within storage of ys
    state = LithBoonkkampIJzermanBracketState(
        b,    # xₙ
        a, # xₙ₋₁
        c,
        fb, # fₙ
        fa, # fₙ₋₁
        fc,
        f′b,
        f′a,
        f′c,
    )

    state
end

function update_state(
    M::LithBoonkkampIJzermanBracket,
    F,
    state::LithBoonkkampIJzermanBracketState{T,S,R},
    options,
    l=NullTracks(),
) where {T,S,R}
    b::T, c::T, a::T = state.xn1, state.c, state.xn0
    fb::S, fc::S, fa::S = state.fxn1, state.fc, state.fxn0
    f′a::R, f′c::R, f′b::R = state.fp0, state.fpc, state.fp1

    # Get next interpolating step
    # decide on S and D;
    # S is 3 if a,b,c are distinct; D=1 unless all derivative info will be of the wrong sign.
    s::Int = ((a == c) || (b == c)) ? 2 : 3

    # which derivatives do we include
    sₘ = sign((fb - fa) / (b - a))
    mc, mb = sign(f′c) == sₘ, sign(f′b) == sₘ

    d₀::S = zero(S)
    if s == 2
        if mc || mb
            # D = 1
            a2s, b2s = lmm_coefficients(LithBoonkkampIJzerman{2,1}(), (c, b), (fc, fb))
            h = -fb

            d₀ = -sum(a2s .* (c, b))
            mb && (d₀ += h * b2s[2] / f′b)
            mc && (d₀ += h * b2s[1] / f′c)
        else
            d₀ = lmm(LithBoonkkampIJzerman{2,0}(), (c, b), (fc, fb))
        end
    else
        ma = sign(f′a) == sₘ
        if mc || mb || ma
            # D = 1
            as, bs = lmm_coefficients(LithBoonkkampIJzerman{3,1}(), (a, c, b), (fa, fc, fb))
            h = -fb

            d₀ = -sum(as .* (a, c, b))
            mb && (d₀ += h * bs[end] / f′b) # only when helpful
            mc && (d₀ += h * bs[end - 1] / f′c)
            ma && (d₀ += h * bs[end - 2] / f′a)
        else
            d₀ = lmm(LithBoonkkampIJzerman{3,0}(), (a, c, b), (fa, fc, fb))
        end
    end

    # If the step is smaller than the tolerance, use the tolerance as step size.
    xatol, xrtol = options.xabstol, options.xreltol
    δ = xatol + abs(b) * xrtol
    Δ₀ = b - d₀
    if abs(Δ₀) <= δ
        d₀ = b - sign(Δ₀) * δ
    end

    # compare to bisection step; extra function evaluation
    d₁ = a + (b - a) * (0.5) #_middle(a, b)
    f₀, Δf₀ = F(d₀)
    f₁, Δf₁ = F(d₁)

    # interpolation outside a,b or bisection better use that
    d::T, fd::S, f′d::S = zero(T), zero(S), zero(S)
    if (abs(f₀) < abs(f₁)) && (min(a, b) < d₀ < max(a, b))
        d, fd, f′d = d₀, f₀, f₀ / Δf₀# interp
    else
        d, fd, f′d = d₁, f₁, f₁ / Δf₁#  bisection
    end

    # either [a,d] a bracket or [d,b]
    # [a < d] < b ...c -- b -> d, c-> b (update?)
    # a < [d < b] ...c -- a -> d (update?)
    if sign(fa) * sign(fd) <= 0
        c, fc, f′c = b, fb, f′b
        b, fb, f′b = d, fd, f′d
    else
        a, fa, f′a = d, fd, f′d
    end

    # a,b bracket; keep |fb| ≤ |fa|
    if abs(fa) < abs(fb)
        c, fc, f′c = b, fb, f′b
        a, b, fa, fb, f′a, f′b = b, a, fb, fa, f′b, f′a
    end

    incfn(l, 3)
    @set! state.xn1 = b
    @set! state.xn0 = a
    @set! state.c = c
    @set! state.fxn1 = fb
    @set! state.fxn0 = fa
    @set! state.fc = fc
    @set! state.fp0 = f′a
    @set! state.fpc = f′c
    @set! state.fp1 = f′b

    return (state, false)
end

function default_tolerances(
    ::M,
    ::Type{T},
    ::Type{S},
) where {M<:LithBoonkkampIJzermanBracket,T,S}
    xatol = 2eps(T)
    xrtol = zero(one(T))
    atol = zero(float(one(S))) * oneunit(S)
    rtol = 2eps(float(one(S))) * one(S)
    maxevals = typemax(Int)
    strict = true
    (xatol, xrtol, atol, rtol, maxevals, strict)
end

### ------
# Script used to generate expressions
#
# Finds expressions by assuming an interpolating polynomial
# goes through the points. From (20) in the paper
# some cases take a **long** time to run.
# At some point, early on, this gets to be more of an academic exercise than
# a performant solution

#= -------
using SymPy

# see Ansatz (20) on p10
function inverse_polynomial_interpretation(s=2,d=2)

    @vars y
    hs = [Sym("h$i") for i ∈ 0:(1+d)*s-1]
    xs = [Sym("x$i") for i ∈ 0:s-1]
    fs = [Sym("f$i") for i ∈ 0:s-1]
    f′s = [Sym("f′$i") for i ∈ 0:s-1]
    f′′s = [Sym("f′′$i") for i ∈ 0:s-1]
    f′′′s = [Sym("f′′′$i") for i ∈ 0:s-1]
    f′′′′s = [Sym("f′′′′$i") for i ∈ 0:s-1]
    f′′′′′s = [Sym("f′′′′′$i") for i ∈ 0:s-1]
    f′′′′′′s = [Sym("f′′′′′′$i") for i ∈ 0:s-1]
    h0 = first(hs)

    H(y) = sum(hs[i+1]*(y-fs[s])^i for i ∈ 0:(1+d)*s-1)
    Hⁱ = H⁰ = H(y)
    Hⁱs = eltype(H⁰)[]
    for i ∈ 1:d
        Hⁱ = diff(Hⁱ,y)
        push!(Hⁱs, Hⁱ)
    end

    eqs = Sym[subs(H(fs[i]), Dict(h0=>xs[s])) - xs[i] for i ∈ 1:s-1]
    for i ∈ 1:s
        # cf. Liptag
        f1,f2,f3,f4,f5,f6 = f′s[i],f′′s[i],f′′′s[i],f′′′′s[i],f′′′′′s[i],f′′′′′′s[i]
        g′ = 1/f1
        g′′ = -f2/f1^3
        g′′′ = (3*f2^2 - f1*f3)/(f1^5)
        g′′′′ = -(15*f2^2 + 10*f1*f2*f3+f1^2*f4)/f1^7
        g′′′′′ = (105*f2^4 -105*f1*f2^2*f3 + 10*f1^2*f3^2 + 15*f1^2*f2*f4 -f1^3*f5)/f1^9
        g′′′′′′ = (-f1^4*f6 + 21*f1^3*f2*f5 + 35*f1^3*f3*f4 - 210*f1^2*f2^2*f4 - 280*f1^2*f2*f3^2 + 1260*f1*f2^3*f3 - 945*f2^5)/f1^11
        gⁱs = [g′,g′′,g′′′, g′′′′,g′′′′′,g′′′′′′]

        for j ∈ 1:d
            push!(eqs, subs(Hⁱs[j], Dict(y=>fs[i], h0=>xs[s])) - gⁱs[j])
        end

    end

    ϕ = sympy.linsolve(eqs, hs[2:end]...)
    ϕ = first(elements(ϕ))
    ϕ = Sym.(convert(Tuple, ϕ.__pyobject__))
    D = Dict{Any,Any}(h0=>xs[s])
    for i in 1:(d+1)*s-1
        D[hs[i+1]] = ϕ[i]
    end
    subs(H(0), D) |> simplify
end

# For g = f⁻¹ return [g', g'', g''',..., g⁽ⁿ⁾]
# (cf [Liptaj](https://vixra.org/pdf/1703.0295v1.pdf)
function liptag(N)

    @vars x₀ Δₓ
    fs = [Sym("f$i") for i ∈ 1:N]
    gs = [Sym("g$i") for i ∈ 1:N]
    a(i) = fs[i]/factorial(i)
    b(i) = gs[i]/factorial(i)

    gᵏs = [1/fs[1]]
    for n ∈ 2:N
        Δy = sum(a(j) * Δₓ^j for j ∈ 1:n)
        l = x₀ + Δₓ
        r = x₀ + sum(b(i)*Δy^i for i ∈ 1:n)
        ϕ = solve(l-r, gs[n])[1]

        for j ∈ 1:n-1
            ϕ = subs(ϕ, gs[j] => gᵏs[j])
        end
        L = limit(ϕ, Δₓ => 0)
        push!(gᵏs, L)
    end
    gᵏs
end
=#

# have computed these
# S/D 0 1 2 3 4  5 6
# 1   x ✓ ✓ ✓ ✓ ✓ ✓
# 2   ✓ ✓ ✓ ✓ ✓ ✓ ✓
# 3   ✓ ✓ - - x x x
# 4   ✓ - x x x x x
# 5   ✓ x x x x x x
# 6   - x x x x x x

# - can be found with script, but answers are too long for
#   inclusion here

## We have two means to do this:
## Using coefficients as,bs, ... returned by lmm_coefficients
## x = ∑ aᵢxᵢ + ∑ⱼ₊₁ⁿ ∑ᵢ bʲᵢFʲᵢ, where Fʲ is the jth derivative of g⁻¹ (F¹ = 1/f'...)
## Using a polynomial interpolant, H(y), going through (xᵢ,fʲ(xᵢ)), j ∈ 0:N)

function lmm(::Val{S}, ::Val{D}, xs, ys) where {S,D}
    xi = ntuple(ii -> xs[ii], Val(S))
    yi = ntuple(ii -> ntuple(j -> ys[ii][j], Val(S)), Val(D + 1))
    lmm(LithBoonkkampIJzerman{S,D}(), xi, yi...)
end

# secant
function lmm(::LithBoonkkampIJzerman{2,0}, xs, fs)
    x0, x1 = xs
    f0, f1 = fs

    (f0 * x1 - f1 * x0) / (f0 - f1)
end

function lmm(::LithBoonkkampIJzerman{3,0}, xs, fs)
    x0, x1, x2 = xs
    f0, f1, f2 = fs

    (
        f0^2 * f1 * x2 - f0^2 * f2 * x1 - f0 * f1^2 * x2 + f0 * f2^2 * x1 + f1^2 * f2 * x0 -
        f1 * f2^2 * x0
    ) / (f0^2 * f1 - f0^2 * f2 - f0 * f1^2 + f0 * f2^2 + f1^2 * f2 - f1 * f2^2)
end

function lmm(::LithBoonkkampIJzerman{4,0}, xs, fs)
    x0, x1, x2, x3 = xs
    f0, f1, f2, f3 = fs

    (
        f0^3 * f1^2 * f2 * x3 - f0^3 * f1^2 * f3 * x2 - f0^3 * f1 * f2^2 * x3 +
        f0^3 * f1 * f3^2 * x2 +
        f0^3 * f2^2 * f3 * x1 - f0^3 * f2 * f3^2 * x1 - f0^2 * f1^3 * f2 * x3 +
        f0^2 * f1^3 * f3 * x2 +
        f0^2 * f1 * f2^3 * x3 - f0^2 * f1 * f3^3 * x2 - f0^2 * f2^3 * f3 * x1 +
        f0^2 * f2 * f3^3 * x1 +
        f0 * f1^3 * f2^2 * x3 - f0 * f1^3 * f3^2 * x2 - f0 * f1^2 * f2^3 * x3 +
        f0 * f1^2 * f3^3 * x2 +
        f0 * f2^3 * f3^2 * x1 - f0 * f2^2 * f3^3 * x1 - f1^3 * f2^2 * f3 * x0 +
        f1^3 * f2 * f3^2 * x0 +
        f1^2 * f2^3 * f3 * x0 - f1^2 * f2 * f3^3 * x0 - f1 * f2^3 * f3^2 * x0 +
        f1 * f2^2 * f3^3 * x0
    ) / (
        f0^3 * f1^2 * f2 - f0^3 * f1^2 * f3 - f0^3 * f1 * f2^2 +
        f0^3 * f1 * f3^2 +
        f0^3 * f2^2 * f3 - f0^3 * f2 * f3^2 - f0^2 * f1^3 * f2 +
        f0^2 * f1^3 * f3 +
        f0^2 * f1 * f2^3 - f0^2 * f1 * f3^3 - f0^2 * f2^3 * f3 +
        f0^2 * f2 * f3^3 +
        f0 * f1^3 * f2^2 - f0 * f1^3 * f3^2 - f0 * f1^2 * f2^3 +
        f0 * f1^2 * f3^3 +
        f0 * f2^3 * f3^2 - f0 * f2^2 * f3^3 - f1^3 * f2^2 * f3 +
        f1^3 * f2 * f3^2 +
        f1^2 * f2^3 * f3 - f1^2 * f2 * f3^3 - f1 * f2^3 * f3^2 + f1 * f2^2 * f3^3
    )
end

function lmm(::LithBoonkkampIJzerman{5,0}, xs, fs)
    x0, x1, x2, x3, x4 = xs
    f0, f1, f2, f3, f4 = fs

    (
        f0^4 * f1^3 * f2^2 * f3 * x4 - f0^4 * f1^3 * f2^2 * f4 * x3 -
        f0^4 * f1^3 * f2 * f3^2 * x4 +
        f0^4 * f1^3 * f2 * f4^2 * x3 +
        f0^4 * f1^3 * f3^2 * f4 * x2 - f0^4 * f1^3 * f3 * f4^2 * x2 -
        f0^4 * f1^2 * f2^3 * f3 * x4 +
        f0^4 * f1^2 * f2^3 * f4 * x3 +
        f0^4 * f1^2 * f2 * f3^3 * x4 - f0^4 * f1^2 * f2 * f4^3 * x3 -
        f0^4 * f1^2 * f3^3 * f4 * x2 +
        f0^4 * f1^2 * f3 * f4^3 * x2 +
        f0^4 * f1 * f2^3 * f3^2 * x4 - f0^4 * f1 * f2^3 * f4^2 * x3 -
        f0^4 * f1 * f2^2 * f3^3 * x4 +
        f0^4 * f1 * f2^2 * f4^3 * x3 +
        f0^4 * f1 * f3^3 * f4^2 * x2 - f0^4 * f1 * f3^2 * f4^3 * x2 -
        f0^4 * f2^3 * f3^2 * f4 * x1 +
        f0^4 * f2^3 * f3 * f4^2 * x1 +
        f0^4 * f2^2 * f3^3 * f4 * x1 - f0^4 * f2^2 * f3 * f4^3 * x1 -
        f0^4 * f2 * f3^3 * f4^2 * x1 + f0^4 * f2 * f3^2 * f4^3 * x1 -
        f0^3 * f1^4 * f2^2 * f3 * x4 +
        f0^3 * f1^4 * f2^2 * f4 * x3 +
        f0^3 * f1^4 * f2 * f3^2 * x4 - f0^3 * f1^4 * f2 * f4^2 * x3 -
        f0^3 * f1^4 * f3^2 * f4 * x2 +
        f0^3 * f1^4 * f3 * f4^2 * x2 +
        f0^3 * f1^2 * f2^4 * f3 * x4 - f0^3 * f1^2 * f2^4 * f4 * x3 -
        f0^3 * f1^2 * f2 * f3^4 * x4 +
        f0^3 * f1^2 * f2 * f4^4 * x3 +
        f0^3 * f1^2 * f3^4 * f4 * x2 - f0^3 * f1^2 * f3 * f4^4 * x2 -
        f0^3 * f1 * f2^4 * f3^2 * x4 +
        f0^3 * f1 * f2^4 * f4^2 * x3 +
        f0^3 * f1 * f2^2 * f3^4 * x4 - f0^3 * f1 * f2^2 * f4^4 * x3 -
        f0^3 * f1 * f3^4 * f4^2 * x2 +
        f0^3 * f1 * f3^2 * f4^4 * x2 +
        f0^3 * f2^4 * f3^2 * f4 * x1 - f0^3 * f2^4 * f3 * f4^2 * x1 -
        f0^3 * f2^2 * f3^4 * f4 * x1 +
        f0^3 * f2^2 * f3 * f4^4 * x1 +
        f0^3 * f2 * f3^4 * f4^2 * x1 - f0^3 * f2 * f3^2 * f4^4 * x1 +
        f0^2 * f1^4 * f2^3 * f3 * x4 - f0^2 * f1^4 * f2^3 * f4 * x3 -
        f0^2 * f1^4 * f2 * f3^3 * x4 +
        f0^2 * f1^4 * f2 * f4^3 * x3 +
        f0^2 * f1^4 * f3^3 * f4 * x2 - f0^2 * f1^4 * f3 * f4^3 * x2 -
        f0^2 * f1^3 * f2^4 * f3 * x4 +
        f0^2 * f1^3 * f2^4 * f4 * x3 +
        f0^2 * f1^3 * f2 * f3^4 * x4 - f0^2 * f1^3 * f2 * f4^4 * x3 -
        f0^2 * f1^3 * f3^4 * f4 * x2 +
        f0^2 * f1^3 * f3 * f4^4 * x2 +
        f0^2 * f1 * f2^4 * f3^3 * x4 - f0^2 * f1 * f2^4 * f4^3 * x3 -
        f0^2 * f1 * f2^3 * f3^4 * x4 +
        f0^2 * f1 * f2^3 * f4^4 * x3 +
        f0^2 * f1 * f3^4 * f4^3 * x2 - f0^2 * f1 * f3^3 * f4^4 * x2 -
        f0^2 * f2^4 * f3^3 * f4 * x1 +
        f0^2 * f2^4 * f3 * f4^3 * x1 +
        f0^2 * f2^3 * f3^4 * f4 * x1 - f0^2 * f2^3 * f3 * f4^4 * x1 -
        f0^2 * f2 * f3^4 * f4^3 * x1 + f0^2 * f2 * f3^3 * f4^4 * x1 -
        f0 * f1^4 * f2^3 * f3^2 * x4 +
        f0 * f1^4 * f2^3 * f4^2 * x3 +
        f0 * f1^4 * f2^2 * f3^3 * x4 - f0 * f1^4 * f2^2 * f4^3 * x3 -
        f0 * f1^4 * f3^3 * f4^2 * x2 +
        f0 * f1^4 * f3^2 * f4^3 * x2 +
        f0 * f1^3 * f2^4 * f3^2 * x4 - f0 * f1^3 * f2^4 * f4^2 * x3 -
        f0 * f1^3 * f2^2 * f3^4 * x4 +
        f0 * f1^3 * f2^2 * f4^4 * x3 +
        f0 * f1^3 * f3^4 * f4^2 * x2 - f0 * f1^3 * f3^2 * f4^4 * x2 -
        f0 * f1^2 * f2^4 * f3^3 * x4 +
        f0 * f1^2 * f2^4 * f4^3 * x3 +
        f0 * f1^2 * f2^3 * f3^4 * x4 - f0 * f1^2 * f2^3 * f4^4 * x3 -
        f0 * f1^2 * f3^4 * f4^3 * x2 +
        f0 * f1^2 * f3^3 * f4^4 * x2 +
        f0 * f2^4 * f3^3 * f4^2 * x1 - f0 * f2^4 * f3^2 * f4^3 * x1 -
        f0 * f2^3 * f3^4 * f4^2 * x1 +
        f0 * f2^3 * f3^2 * f4^4 * x1 +
        f0 * f2^2 * f3^4 * f4^3 * x1 - f0 * f2^2 * f3^3 * f4^4 * x1 +
        f1^4 * f2^3 * f3^2 * f4 * x0 - f1^4 * f2^3 * f3 * f4^2 * x0 -
        f1^4 * f2^2 * f3^3 * f4 * x0 +
        f1^4 * f2^2 * f3 * f4^3 * x0 +
        f1^4 * f2 * f3^3 * f4^2 * x0 - f1^4 * f2 * f3^2 * f4^3 * x0 -
        f1^3 * f2^4 * f3^2 * f4 * x0 +
        f1^3 * f2^4 * f3 * f4^2 * x0 +
        f1^3 * f2^2 * f3^4 * f4 * x0 - f1^3 * f2^2 * f3 * f4^4 * x0 -
        f1^3 * f2 * f3^4 * f4^2 * x0 +
        f1^3 * f2 * f3^2 * f4^4 * x0 +
        f1^2 * f2^4 * f3^3 * f4 * x0 - f1^2 * f2^4 * f3 * f4^3 * x0 -
        f1^2 * f2^3 * f3^4 * f4 * x0 +
        f1^2 * f2^3 * f3 * f4^4 * x0 +
        f1^2 * f2 * f3^4 * f4^3 * x0 - f1^2 * f2 * f3^3 * f4^4 * x0 -
        f1 * f2^4 * f3^3 * f4^2 * x0 +
        f1 * f2^4 * f3^2 * f4^3 * x0 +
        f1 * f2^3 * f3^4 * f4^2 * x0 - f1 * f2^3 * f3^2 * f4^4 * x0 -
        f1 * f2^2 * f3^4 * f4^3 * x0 + f1 * f2^2 * f3^3 * f4^4 * x0
    ) / (
        f0^4 * f1^3 * f2^2 * f3 - f0^4 * f1^3 * f2^2 * f4 - f0^4 * f1^3 * f2 * f3^2 +
        f0^4 * f1^3 * f2 * f4^2 +
        f0^4 * f1^3 * f3^2 * f4 - f0^4 * f1^3 * f3 * f4^2 - f0^4 * f1^2 * f2^3 * f3 +
        f0^4 * f1^2 * f2^3 * f4 +
        f0^4 * f1^2 * f2 * f3^3 - f0^4 * f1^2 * f2 * f4^3 - f0^4 * f1^2 * f3^3 * f4 +
        f0^4 * f1^2 * f3 * f4^3 +
        f0^4 * f1 * f2^3 * f3^2 - f0^4 * f1 * f2^3 * f4^2 - f0^4 * f1 * f2^2 * f3^3 +
        f0^4 * f1 * f2^2 * f4^3 +
        f0^4 * f1 * f3^3 * f4^2 - f0^4 * f1 * f3^2 * f4^3 - f0^4 * f2^3 * f3^2 * f4 +
        f0^4 * f2^3 * f3 * f4^2 +
        f0^4 * f2^2 * f3^3 * f4 - f0^4 * f2^2 * f3 * f4^3 - f0^4 * f2 * f3^3 * f4^2 +
        f0^4 * f2 * f3^2 * f4^3 - f0^3 * f1^4 * f2^2 * f3 +
        f0^3 * f1^4 * f2^2 * f4 +
        f0^3 * f1^4 * f2 * f3^2 - f0^3 * f1^4 * f2 * f4^2 - f0^3 * f1^4 * f3^2 * f4 +
        f0^3 * f1^4 * f3 * f4^2 +
        f0^3 * f1^2 * f2^4 * f3 - f0^3 * f1^2 * f2^4 * f4 - f0^3 * f1^2 * f2 * f3^4 +
        f0^3 * f1^2 * f2 * f4^4 +
        f0^3 * f1^2 * f3^4 * f4 - f0^3 * f1^2 * f3 * f4^4 - f0^3 * f1 * f2^4 * f3^2 +
        f0^3 * f1 * f2^4 * f4^2 +
        f0^3 * f1 * f2^2 * f3^4 - f0^3 * f1 * f2^2 * f4^4 - f0^3 * f1 * f3^4 * f4^2 +
        f0^3 * f1 * f3^2 * f4^4 +
        f0^3 * f2^4 * f3^2 * f4 - f0^3 * f2^4 * f3 * f4^2 - f0^3 * f2^2 * f3^4 * f4 +
        f0^3 * f2^2 * f3 * f4^4 +
        f0^3 * f2 * f3^4 * f4^2 - f0^3 * f2 * f3^2 * f4^4 + f0^2 * f1^4 * f2^3 * f3 -
        f0^2 * f1^4 * f2^3 * f4 - f0^2 * f1^4 * f2 * f3^3 +
        f0^2 * f1^4 * f2 * f4^3 +
        f0^2 * f1^4 * f3^3 * f4 - f0^2 * f1^4 * f3 * f4^3 - f0^2 * f1^3 * f2^4 * f3 +
        f0^2 * f1^3 * f2^4 * f4 +
        f0^2 * f1^3 * f2 * f3^4 - f0^2 * f1^3 * f2 * f4^4 - f0^2 * f1^3 * f3^4 * f4 +
        f0^2 * f1^3 * f3 * f4^4 +
        f0^2 * f1 * f2^4 * f3^3 - f0^2 * f1 * f2^4 * f4^3 - f0^2 * f1 * f2^3 * f3^4 +
        f0^2 * f1 * f2^3 * f4^4 +
        f0^2 * f1 * f3^4 * f4^3 - f0^2 * f1 * f3^3 * f4^4 - f0^2 * f2^4 * f3^3 * f4 +
        f0^2 * f2^4 * f3 * f4^3 +
        f0^2 * f2^3 * f3^4 * f4 - f0^2 * f2^3 * f3 * f4^4 - f0^2 * f2 * f3^4 * f4^3 +
        f0^2 * f2 * f3^3 * f4^4 - f0 * f1^4 * f2^3 * f3^2 +
        f0 * f1^4 * f2^3 * f4^2 +
        f0 * f1^4 * f2^2 * f3^3 - f0 * f1^4 * f2^2 * f4^3 - f0 * f1^4 * f3^3 * f4^2 +
        f0 * f1^4 * f3^2 * f4^3 +
        f0 * f1^3 * f2^4 * f3^2 - f0 * f1^3 * f2^4 * f4^2 - f0 * f1^3 * f2^2 * f3^4 +
        f0 * f1^3 * f2^2 * f4^4 +
        f0 * f1^3 * f3^4 * f4^2 - f0 * f1^3 * f3^2 * f4^4 - f0 * f1^2 * f2^4 * f3^3 +
        f0 * f1^2 * f2^4 * f4^3 +
        f0 * f1^2 * f2^3 * f3^4 - f0 * f1^2 * f2^3 * f4^4 - f0 * f1^2 * f3^4 * f4^3 +
        f0 * f1^2 * f3^3 * f4^4 +
        f0 * f2^4 * f3^3 * f4^2 - f0 * f2^4 * f3^2 * f4^3 - f0 * f2^3 * f3^4 * f4^2 +
        f0 * f2^3 * f3^2 * f4^4 +
        f0 * f2^2 * f3^4 * f4^3 - f0 * f2^2 * f3^3 * f4^4 + f1^4 * f2^3 * f3^2 * f4 -
        f1^4 * f2^3 * f3 * f4^2 - f1^4 * f2^2 * f3^3 * f4 +
        f1^4 * f2^2 * f3 * f4^3 +
        f1^4 * f2 * f3^3 * f4^2 - f1^4 * f2 * f3^2 * f4^3 - f1^3 * f2^4 * f3^2 * f4 +
        f1^3 * f2^4 * f3 * f4^2 +
        f1^3 * f2^2 * f3^4 * f4 - f1^3 * f2^2 * f3 * f4^4 - f1^3 * f2 * f3^4 * f4^2 +
        f1^3 * f2 * f3^2 * f4^4 +
        f1^2 * f2^4 * f3^3 * f4 - f1^2 * f2^4 * f3 * f4^3 - f1^2 * f2^3 * f3^4 * f4 +
        f1^2 * f2^3 * f3 * f4^4 +
        f1^2 * f2 * f3^4 * f4^3 - f1^2 * f2 * f3^3 * f4^4 - f1 * f2^4 * f3^3 * f4^2 +
        f1 * f2^4 * f3^2 * f4^3 +
        f1 * f2^3 * f3^4 * f4^2 - f1 * f2^3 * f3^2 * f4^4 - f1 * f2^2 * f3^4 * f4^3 +
        f1 * f2^2 * f3^3 * f4^4
    )
end

function lmm(::LithBoonkkampIJzerman{6,0}, xs, fs)
    x0, x1, x2, x3, x4, x5 = xs
    f0, f1, f2, f3, f4, f5 = fs

    error("not implemented")
end

## d = 1; Newton-like

# return (as, bs⁰,[bs¹,...,bsⁿ⁻¹])
# where -∑ aᵢ xᵢ + h ⋅ ∑ₙ (∑ bsʲᵢ Fʲᵢ)l
function lmm_coefficients(::LithBoonkkampIJzerman{1,1}, xs, fs)
    a0 = -one(xs[1])
    b0 = one(fs[1])

    return (a0,), (b0,)
end

function lmm_coefficients(::LithBoonkkampIJzerman{2,1}, xs, fs)
    q = fs[1] / fs[2]

    # from the paper
    # x2 + a1 x1 + a0x0 =  h3 * (b1 * 1/fp1 + b0 * 1/fp0)
    a0 = (1 - 3q) / (q - 1)^3
    a1 = -1 - a0
    b0 = q / (q - 1)^2
    b1 = q * b0

    return (a0, a1), (b0, b1)
end

function lmm_coefficients(::LithBoonkkampIJzerman{3,1}, xs, fs)

    # from the paper
    q0 = fs[3 - 2] / fs[3]
    q1 = fs[3 - 1] / fs[3]

    a0 = (q1^2 * (q0 * (3 + 3q1 - 5q0) - q1)) / ((q0 - 1)^3 * (q0 - q1)^3)
    a1 = (q0^2 * (q1 * (5q1 - 3q0 - 3) + q0)) / ((q1 - 1)^3 * (q0 - q1)^3)
    a2 = (q0^2 * q1^2 * (3q1 - q0 * (q1 - 3) - 5)) / ((q0 - 1)^3 * (q1 - 1)^3) # minor typo in (27c)

    b0 = (q0 * q1^2) / ((q0 - 1)^2 * (q0 - q1)^2)
    b1 = (q0^2 * q1) / ((q0 - q1)^2 * (q1 - 1)^2)
    b2 = (q0^2 * q1^2) / ((q0 - 1)^2 * (q1 - 1)^2)

    return (a0, a1, a2), (b0, b1, b2)
end

function lmm_coefficients(::LithBoonkkampIJzerman{S,1}, xs, fs) where {S}
    error("not computed")
end

function lmm(L::LithBoonkkampIJzerman{S,1}, xs, fs, f′s) where {S}
    as, bs = lmm_coefficients(L, xs, fs)
    Fs = 1 ./ f′s # F = (g⁻¹)'
    h = -fs[S]

    -sum(as[i] * xs[i] for i in 1:S) + h * sum(bs[i] * Fs[i] for i in 1:S)
end

function lmm(::LithBoonkkampIJzerman{4,1}, xs, fs, f′s)
    x0, x1, x2, x3 = xs
    f0, f1, f2, f3 = fs
    f′0, f′1, f′2, f′3 = f′s

    # can get with script, but too long as found
    error("not implemented")
end

## d = 2; Halley-like
function lmm(::LithBoonkkampIJzerman{1,2}, xs, fs, f′s, f′′s)
    x0 = xs[1]
    f0 = fs[1]
    f′0 = f′s[1]
    f′′0 = f′′s[1]

    -f0^2 * f′′0 / (2 * f′0^3) - f0 / f′0 + x0
end

function lmm(::LithBoonkkampIJzerman{2,2}, xs, fs, f′s, f′′s)
    x0, x1 = xs
    f0, f1 = fs
    f′0, f′1 = f′s
    f′′0, f′′1 = f′′s

    (
        -f0^5 * f1^2 * f′0^3 * f′′1 / 2 - f0^5 * f1 * f′0^3 * f′1^2 +
        f0^5 * f′0^3 * f′1^3 * x1 +
        f0^4 * f1^3 * f′0^3 * f′′1 +
        f0^4 * f1^3 * f′1^3 * f′′0 / 2 +
        5 * f0^4 * f1^2 * f′0^3 * f′1^2 - 5 * f0^4 * f1 * f′0^3 * f′1^3 * x1 -
        f0^3 * f1^4 * f′0^3 * f′′1 / 2 - f0^3 * f1^4 * f′1^3 * f′′0 -
        4 * f0^3 * f1^3 * f′0^3 * f′1^2 +
        4 * f0^3 * f1^3 * f′0^2 * f′1^3 +
        10 * f0^3 * f1^2 * f′0^3 * f′1^3 * x1 +
        f0^2 * f1^5 * f′1^3 * f′′0 / 2 - 5 * f0^2 * f1^4 * f′0^2 * f′1^3 -
        10 * f0^2 * f1^3 * f′0^3 * f′1^3 * x0 +
        f0 * f1^5 * f′0^2 * f′1^3 +
        5 * f0 * f1^4 * f′0^3 * f′1^3 * x0 - f1^5 * f′0^3 * f′1^3 * x0
    ) / (
        f′0^3 *
        f′1^3 *
        (f0^5 - 5 * f0^4 * f1 + 10 * f0^3 * f1^2 - 10 * f0^2 * f1^3 + 5 * f0 * f1^4 - f1^5)
    )
end

function lmm(::LithBoonkkampIJzerman{3,2}, xs, fs, f′s, f′′s)
    x0, x1, x2 = xs
    f0, f1, f2 = fs
    f′0, f′1, f′2 = f′s
    f′′0, f′′1, f′′2 = f′′s

    ## can get from script, but too long for inclusion here
    error("not implemented")
end

function lmm(::LithBoonkkampIJzerman{4,2}, xs, fs, f′s, f′′s)
    error("Not computed")
end

## d = 3
function lmm(::LithBoonkkampIJzerman{1,3}, xs, fs, f′s, f′′s, f′′′s)
    x0 = xs[1]
    f0 = fs[1]
    f′0 = f′s[1]
    f′′0 = f′′s[1]
    f′′′0 = f′′′s[1]

    (
        f0^3 * (f′0 * f′′′0 - 3 * f′′0^2) / 6 - f0^2 * f′0^2 * f′′0 / 2 - f0 * f′0^4 +
        f′0^5 * x0
    ) / f′0^5
end

function lmm(::LithBoonkkampIJzerman{2,3}, xs, fs, f′s, f′′s, f′′′s)
    x0, x1 = xs
    f0, f1 = fs
    f′0, f′1 = f′s
    f′′0, f′′1 = f′′s
    f′′′0, f′′′1 = f′′′s

    (
        f0^7 * f1^3 * f′0^5 * f′1 * f′′′1 - 3 * f0^7 * f1^3 * f′0^5 * f′′1^2 -
        3 * f0^7 * f1^2 * f′0^5 * f′1^2 * f′′1 - 6 * f0^7 * f1 * f′0^5 * f′1^4 +
        6 * f0^7 * f′0^5 * f′1^5 * x1 - 3 * f0^6 * f1^4 * f′0^5 * f′1 * f′′′1 +
        9 * f0^6 * f1^4 * f′0^5 * f′′1^2 +
        f0^6 * f1^4 * f′0 * f′1^5 * f′′′0 - 3 * f0^6 * f1^4 * f′1^5 * f′′0^2 +
        21 * f0^6 * f1^3 * f′0^5 * f′1^2 * f′′1 +
        42 * f0^6 * f1^2 * f′0^5 * f′1^4 - 42 * f0^6 * f1 * f′0^5 * f′1^5 * x1 +
        3 * f0^5 * f1^5 * f′0^5 * f′1 * f′′′1 - 9 * f0^5 * f1^5 * f′0^5 * f′′1^2 -
        3 * f0^5 * f1^5 * f′0 * f′1^5 * f′′′0 + 9 * f0^5 * f1^5 * f′1^5 * f′′0^2 -
        33 * f0^5 * f1^4 * f′0^5 * f′1^2 * f′′1 - 15 * f0^5 * f1^4 * f′0^2 * f′1^5 * f′′0 -
        126 * f0^5 * f1^3 * f′0^5 * f′1^4 + 126 * f0^5 * f1^2 * f′0^5 * f′1^5 * x1 -
        f0^4 * f1^6 * f′0^5 * f′1 * f′′′1 +
        3 * f0^4 * f1^6 * f′0^5 * f′′1^2 +
        3 * f0^4 * f1^6 * f′0 * f′1^5 * f′′′0 - 9 * f0^4 * f1^6 * f′1^5 * f′′0^2 +
        15 * f0^4 * f1^5 * f′0^5 * f′1^2 * f′′1 +
        33 * f0^4 * f1^5 * f′0^2 * f′1^5 * f′′0 +
        90 * f0^4 * f1^4 * f′0^5 * f′1^4 - 90 * f0^4 * f1^4 * f′0^4 * f′1^5 -
        210 * f0^4 * f1^3 * f′0^5 * f′1^5 * x1 - f0^3 * f1^7 * f′0 * f′1^5 * f′′′0 +
        3 * f0^3 * f1^7 * f′1^5 * f′′0^2 - 21 * f0^3 * f1^6 * f′0^2 * f′1^5 * f′′0 +
        126 * f0^3 * f1^5 * f′0^4 * f′1^5 +
        210 * f0^3 * f1^4 * f′0^5 * f′1^5 * x0 +
        3 * f0^2 * f1^7 * f′0^2 * f′1^5 * f′′0 - 42 * f0^2 * f1^6 * f′0^4 * f′1^5 -
        126 * f0^2 * f1^5 * f′0^5 * f′1^5 * x0 +
        6 * f0 * f1^7 * f′0^4 * f′1^5 +
        42 * f0 * f1^6 * f′0^5 * f′1^5 * x0 - 6 * f1^7 * f′0^5 * f′1^5 * x0
    ) / (
        6 *
        f′0^5 *
        f′1^5 *
        (
            f0^7 - 7 * f0^6 * f1 + 21 * f0^5 * f1^2 - 35 * f0^4 * f1^3 + 35 * f0^3 * f1^4 -
            21 * f0^2 * f1^5 + 7 * f0 * f1^6 - f1^7
        )
    )
end

function lmm(::LithBoonkkampIJzerman{3,3}, xs, fs, f′s, f′′s, f′′′s)
    x0, x1, x2 = xs
    f0, f1, f2 = fs
    f′0, f′1, f′2 = f′s
    f′′0, f′′1, f′′2 = f′′s
    f′′′0, f′′′1, f′′′2 = f′′′s

    # can get from script, but too long for inclusion here
    error("not implemented")
end

function lmm(::LithBoonkkampIJzerman{4,3}, xs, fs, f′s, f′′s, f′′′s)
    x0, x1, x2, x3 = xs
    f0, f1, f2, f3 = fs
    f′0, f′1, f′2, f′3 = f′s
    f′′0, f′′1, f′′2, f′′3 = f′′s
    f′′′0, f′′′1, f′′′2, f′′′3 = f′′′s

    error("not computed")
end

## d = 4
function lmm(::LithBoonkkampIJzerman{1,4}, xs, fs, f′s, f′′s, f′′′s, f′′′′s)
    x0 = xs[1]
    f0 = fs[1]
    f′0 = f′s[1]
    f′′0 = f′′s[1]
    f′′′0 = f′′′s[1]
    f′′′′0 = f′′′′s[1]

    (
        -f0^4 * (f′0^2 * f′′′′0 + 10 * f′0 * f′′0 * f′′′0 + 15 * f′′0^2) / 24 +
        f0^3 * f′0^2 * (f′0 * f′′′0 - 3 * f′′0^2) / 6 - f0^2 * f′0^4 * f′′0 / 2 -
        f0 * f′0^6 + f′0^7 * x0
    ) / f′0^7
end

function lmm(::LithBoonkkampIJzerman{2,4}, xs, fs, f′s, f′′s, f′′′s, f′′′′s)
    x0, x1 = xs
    f0, f1 = fs
    f′0, f′1 = f′s
    f′′0, f′′1 = f′′s
    f′′′0, f′′′1 = f′′′s
    f′′′′0, f′′′′1 = f′′′′s

    (
        -f0^9 * f1^4 * f′0^7 * f′1^2 * f′′′′1 -
        10 * f0^9 * f1^4 * f′0^7 * f′1 * f′′1 * f′′′1 - 15 * f0^9 * f1^4 * f′0^7 * f′′1^2 +
        4 * f0^9 * f1^3 * f′0^7 * f′1^3 * f′′′1 -
        12 * f0^9 * f1^3 * f′0^7 * f′1^2 * f′′1^2 -
        12 * f0^9 * f1^2 * f′0^7 * f′1^4 * f′′1 - 24 * f0^9 * f1 * f′0^7 * f′1^6 +
        24 * f0^9 * f′0^7 * f′1^7 * x1 +
        4 * f0^8 * f1^5 * f′0^7 * f′1^2 * f′′′′1 +
        40 * f0^8 * f1^5 * f′0^7 * f′1 * f′′1 * f′′′1 +
        60 * f0^8 * f1^5 * f′0^7 * f′′1^2 +
        f0^8 * f1^5 * f′0^2 * f′1^7 * f′′′′0 +
        10 * f0^8 * f1^5 * f′0 * f′1^7 * f′′0 * f′′′0 +
        15 * f0^8 * f1^5 * f′1^7 * f′′0^2 - 36 * f0^8 * f1^4 * f′0^7 * f′1^3 * f′′′1 +
        108 * f0^8 * f1^4 * f′0^7 * f′1^2 * f′′1^2 +
        108 * f0^8 * f1^3 * f′0^7 * f′1^4 * f′′1 +
        216 * f0^8 * f1^2 * f′0^7 * f′1^6 - 216 * f0^8 * f1 * f′0^7 * f′1^7 * x1 -
        6 * f0^7 * f1^6 * f′0^7 * f′1^2 * f′′′′1 -
        60 * f0^7 * f1^6 * f′0^7 * f′1 * f′′1 * f′′′1 - 90 * f0^7 * f1^6 * f′0^7 * f′′1^2 -
        4 * f0^7 * f1^6 * f′0^2 * f′1^7 * f′′′′0 -
        40 * f0^7 * f1^6 * f′0 * f′1^7 * f′′0 * f′′′0 - 60 * f0^7 * f1^6 * f′1^7 * f′′0^2 +
        84 * f0^7 * f1^5 * f′0^7 * f′1^3 * f′′′1 -
        252 * f0^7 * f1^5 * f′0^7 * f′1^2 * f′′1^2 -
        24 * f0^7 * f1^5 * f′0^3 * f′1^7 * f′′′0 +
        72 * f0^7 * f1^5 * f′0^2 * f′1^7 * f′′0^2 -
        432 * f0^7 * f1^4 * f′0^7 * f′1^4 * f′′1 - 864 * f0^7 * f1^3 * f′0^7 * f′1^6 +
        864 * f0^7 * f1^2 * f′0^7 * f′1^7 * x1 +
        4 * f0^6 * f1^7 * f′0^7 * f′1^2 * f′′′′1 +
        40 * f0^6 * f1^7 * f′0^7 * f′1 * f′′1 * f′′′1 +
        60 * f0^6 * f1^7 * f′0^7 * f′′1^2 +
        6 * f0^6 * f1^7 * f′0^2 * f′1^7 * f′′′′0 +
        60 * f0^6 * f1^7 * f′0 * f′1^7 * f′′0 * f′′′0 +
        90 * f0^6 * f1^7 * f′1^7 * f′′0^2 - 76 * f0^6 * f1^6 * f′0^7 * f′1^3 * f′′′1 +
        228 * f0^6 * f1^6 * f′0^7 * f′1^2 * f′′1^2 +
        76 * f0^6 * f1^6 * f′0^3 * f′1^7 * f′′′0 -
        228 * f0^6 * f1^6 * f′0^2 * f′1^7 * f′′0^2 +
        588 * f0^6 * f1^5 * f′0^7 * f′1^4 * f′′1 +
        252 * f0^6 * f1^5 * f′0^4 * f′1^7 * f′′0 +
        2016 * f0^6 * f1^4 * f′0^7 * f′1^6 - 2016 * f0^6 * f1^3 * f′0^7 * f′1^7 * x1 -
        f0^5 * f1^8 * f′0^7 * f′1^2 * f′′′′1 -
        10 * f0^5 * f1^8 * f′0^7 * f′1 * f′′1 * f′′′1 - 15 * f0^5 * f1^8 * f′0^7 * f′′1^2 -
        4 * f0^5 * f1^8 * f′0^2 * f′1^7 * f′′′′0 -
        40 * f0^5 * f1^8 * f′0 * f′1^7 * f′′0 * f′′′0 - 60 * f0^5 * f1^8 * f′1^7 * f′′0^2 +
        24 * f0^5 * f1^7 * f′0^7 * f′1^3 * f′′′1 -
        72 * f0^5 * f1^7 * f′0^7 * f′1^2 * f′′1^2 -
        84 * f0^5 * f1^7 * f′0^3 * f′1^7 * f′′′0 +
        252 * f0^5 * f1^7 * f′0^2 * f′1^7 * f′′0^2 -
        252 * f0^5 * f1^6 * f′0^7 * f′1^4 * f′′1 -
        588 * f0^5 * f1^6 * f′0^4 * f′1^7 * f′′0 - 1344 * f0^5 * f1^5 * f′0^7 * f′1^6 +
        1344 * f0^5 * f1^5 * f′0^6 * f′1^7 +
        3024 * f0^5 * f1^4 * f′0^7 * f′1^7 * x1 +
        f0^4 * f1^9 * f′0^2 * f′1^7 * f′′′′0 +
        10 * f0^4 * f1^9 * f′0 * f′1^7 * f′′0 * f′′′0 +
        15 * f0^4 * f1^9 * f′1^7 * f′′0^2 +
        36 * f0^4 * f1^8 * f′0^3 * f′1^7 * f′′′0 -
        108 * f0^4 * f1^8 * f′0^2 * f′1^7 * f′′0^2 +
        432 * f0^4 * f1^7 * f′0^4 * f′1^7 * f′′0 - 2016 * f0^4 * f1^6 * f′0^6 * f′1^7 -
        3024 * f0^4 * f1^5 * f′0^7 * f′1^7 * x0 - 4 * f0^3 * f1^9 * f′0^3 * f′1^7 * f′′′0 +
        12 * f0^3 * f1^9 * f′0^2 * f′1^7 * f′′0^2 -
        108 * f0^3 * f1^8 * f′0^4 * f′1^7 * f′′0 +
        864 * f0^3 * f1^7 * f′0^6 * f′1^7 +
        2016 * f0^3 * f1^6 * f′0^7 * f′1^7 * x0 +
        12 * f0^2 * f1^9 * f′0^4 * f′1^7 * f′′0 - 216 * f0^2 * f1^8 * f′0^6 * f′1^7 -
        864 * f0^2 * f1^7 * f′0^7 * f′1^7 * x0 +
        24 * f0 * f1^9 * f′0^6 * f′1^7 +
        216 * f0 * f1^8 * f′0^7 * f′1^7 * x0 - 24 * f1^9 * f′0^7 * f′1^7 * x0
    ) / (
        24 *
        f′0^7 *
        f′1^7 *
        (
            f0^9 - 9 * f0^8 * f1 + 36 * f0^7 * f1^2 - 84 * f0^6 * f1^3 + 126 * f0^5 * f1^4 -
            126 * f0^4 * f1^5 + 84 * f0^3 * f1^6 - 36 * f0^2 * f1^7 + 9 * f0 * f1^8 - f1^9
        )
    )
end

function lmm(::LithBoonkkampIJzerman{3,4}, xs, fs, f′s, f′′s, f′′′s, f′′′′s)
    x0, x1, x2 = xs
    f0, f1, f2 = fs
    f′0, f′1, f′2 = f′s
    f′′0, f′′1, f′′2 = f′′s
    f′′′0, f′′′1, f′′′2 = f′′′s
    f′′′′0, f′′′′1, f′′′′2 = f′′′′s

    error("not computed")
end

# n = 5

function lmm(::LithBoonkkampIJzerman{1,5}, xs, fs, f′s, f′′s, f′′′s, f′′′′s, f′′′′′s)
    x0 = xs[1]
    f0 = fs[1]
    f′0 = f′s[1]
    f′′0 = f′′s[1]
    f′′′0 = f′′′s[1]
    f′′′′0 = f′′′′s[1]
    f′′′′′0 = f′′′′′s[1]

    (
        f0^5 * (
            f′0^3 * f′′′′′0 - 15 * f′0^2 * f′′0 * f′′′′0 - 10 * f′0^2 * f′′′0^2 +
            105 * f′0 * f′′0^2 * f′′′0 - 105 * f′′0^4
        ) - 5 * f0^4 * f′0^2 * (f′0^2 * f′′′′0 + 10 * f′0 * f′′0 * f′′′0 + 15 * f′′0^2) +
        20 * f0^3 * f′0^4 * (f′0 * f′′′0 - 3 * f′′0^2) - 60 * f0^2 * f′0^6 * f′′0 -
        120 * f0 * f′0^8 + 120 * f′0^9 * x0
    ) / (120 * f′0^9)
end

function lmm(::LithBoonkkampIJzerman{2,5}, xs, fs, f′s, f′′s, f′′′s, f′′′′s, f′′′′′s)
    x0, x1 = xs
    f0, f1 = fs
    f′0, f′1 = f′s
    f′′0, f′′1 = f′′s
    f′′′0, f′′′1 = f′′′s
    f′′′′0, f′′′′1 = f′′′′s
    f′′′′′0, f′′′′′1 = f′′′′′s

    (
        f0^11 * f1^5 * f′0^9 * f′1^3 * f′′′′′1 -
        15 * f0^11 * f1^5 * f′0^9 * f′1^2 * f′′1 * f′′′′1 -
        10 * f0^11 * f1^5 * f′0^9 * f′1^2 * f′′′1^2 +
        105 * f0^11 * f1^5 * f′0^9 * f′1 * f′′1^2 * f′′′1 -
        105 * f0^11 * f1^5 * f′0^9 * f′′1^4 - 5 * f0^11 * f1^4 * f′0^9 * f′1^4 * f′′′′1 -
        50 * f0^11 * f1^4 * f′0^9 * f′1^3 * f′′1 * f′′′1 -
        75 * f0^11 * f1^4 * f′0^9 * f′1^2 * f′′1^2 +
        20 * f0^11 * f1^3 * f′0^9 * f′1^5 * f′′′1 -
        60 * f0^11 * f1^3 * f′0^9 * f′1^4 * f′′1^2 -
        60 * f0^11 * f1^2 * f′0^9 * f′1^6 * f′′1 - 120 * f0^11 * f1 * f′0^9 * f′1^8 +
        120 * f0^11 * f′0^9 * f′1^9 * x1 - 5 * f0^10 * f1^6 * f′0^9 * f′1^3 * f′′′′′1 +
        75 * f0^10 * f1^6 * f′0^9 * f′1^2 * f′′1 * f′′′′1 +
        50 * f0^10 * f1^6 * f′0^9 * f′1^2 * f′′′1^2 -
        525 * f0^10 * f1^6 * f′0^9 * f′1 * f′′1^2 * f′′′1 +
        525 * f0^10 * f1^6 * f′0^9 * f′′1^4 +
        f0^10 * f1^6 * f′0^3 * f′1^9 * f′′′′′0 -
        15 * f0^10 * f1^6 * f′0^2 * f′1^9 * f′′0 * f′′′′0 -
        10 * f0^10 * f1^6 * f′0^2 * f′1^9 * f′′′0^2 +
        105 * f0^10 * f1^6 * f′0 * f′1^9 * f′′0^2 * f′′′0 -
        105 * f0^10 * f1^6 * f′1^9 * f′′0^4 +
        55 * f0^10 * f1^5 * f′0^9 * f′1^4 * f′′′′1 +
        550 * f0^10 * f1^5 * f′0^9 * f′1^3 * f′′1 * f′′′1 +
        825 * f0^10 * f1^5 * f′0^9 * f′1^2 * f′′1^2 -
        220 * f0^10 * f1^4 * f′0^9 * f′1^5 * f′′′1 +
        660 * f0^10 * f1^4 * f′0^9 * f′1^4 * f′′1^2 +
        660 * f0^10 * f1^3 * f′0^9 * f′1^6 * f′′1 +
        1320 * f0^10 * f1^2 * f′0^9 * f′1^8 - 1320 * f0^10 * f1 * f′0^9 * f′1^9 * x1 +
        10 * f0^9 * f1^7 * f′0^9 * f′1^3 * f′′′′′1 -
        150 * f0^9 * f1^7 * f′0^9 * f′1^2 * f′′1 * f′′′′1 -
        100 * f0^9 * f1^7 * f′0^9 * f′1^2 * f′′′1^2 +
        1050 * f0^9 * f1^7 * f′0^9 * f′1 * f′′1^2 * f′′′1 -
        1050 * f0^9 * f1^7 * f′0^9 * f′′1^4 - 5 * f0^9 * f1^7 * f′0^3 * f′1^9 * f′′′′′0 +
        75 * f0^9 * f1^7 * f′0^2 * f′1^9 * f′′0 * f′′′′0 +
        50 * f0^9 * f1^7 * f′0^2 * f′1^9 * f′′′0^2 -
        525 * f0^9 * f1^7 * f′0 * f′1^9 * f′′0^2 * f′′′0 +
        525 * f0^9 * f1^7 * f′1^9 * f′′0^4 - 170 * f0^9 * f1^6 * f′0^9 * f′1^4 * f′′′′1 -
        1700 * f0^9 * f1^6 * f′0^9 * f′1^3 * f′′1 * f′′′1 -
        2550 * f0^9 * f1^6 * f′0^9 * f′1^2 * f′′1^2 -
        35 * f0^9 * f1^6 * f′0^4 * f′1^9 * f′′′′0 -
        350 * f0^9 * f1^6 * f′0^3 * f′1^9 * f′′0 * f′′′0 -
        525 * f0^9 * f1^6 * f′0^2 * f′1^9 * f′′0^2 +
        1100 * f0^9 * f1^5 * f′0^9 * f′1^5 * f′′′1 -
        3300 * f0^9 * f1^5 * f′0^9 * f′1^4 * f′′1^2 -
        3300 * f0^9 * f1^4 * f′0^9 * f′1^6 * f′′1 - 6600 * f0^9 * f1^3 * f′0^9 * f′1^8 +
        6600 * f0^9 * f1^2 * f′0^9 * f′1^9 * x1 -
        10 * f0^8 * f1^8 * f′0^9 * f′1^3 * f′′′′′1 +
        150 * f0^8 * f1^8 * f′0^9 * f′1^2 * f′′1 * f′′′′1 +
        100 * f0^8 * f1^8 * f′0^9 * f′1^2 * f′′′1^2 -
        1050 * f0^8 * f1^8 * f′0^9 * f′1 * f′′1^2 * f′′′1 +
        1050 * f0^8 * f1^8 * f′0^9 * f′′1^4 +
        10 * f0^8 * f1^8 * f′0^3 * f′1^9 * f′′′′′0 -
        150 * f0^8 * f1^8 * f′0^2 * f′1^9 * f′′0 * f′′′′0 -
        100 * f0^8 * f1^8 * f′0^2 * f′1^9 * f′′′0^2 +
        1050 * f0^8 * f1^8 * f′0 * f′1^9 * f′′0^2 * f′′′0 -
        1050 * f0^8 * f1^8 * f′1^9 * f′′0^4 +
        230 * f0^8 * f1^7 * f′0^9 * f′1^4 * f′′′′1 +
        2300 * f0^8 * f1^7 * f′0^9 * f′1^3 * f′′1 * f′′′1 +
        3450 * f0^8 * f1^7 * f′0^9 * f′1^2 * f′′1^2 +
        145 * f0^8 * f1^7 * f′0^4 * f′1^9 * f′′′′0 +
        1450 * f0^8 * f1^7 * f′0^3 * f′1^9 * f′′0 * f′′′0 +
        2175 * f0^8 * f1^7 * f′0^2 * f′1^9 * f′′0^2 -
        2180 * f0^8 * f1^6 * f′0^9 * f′1^5 * f′′′1 +
        6540 * f0^8 * f1^6 * f′0^9 * f′1^4 * f′′1^2 +
        560 * f0^8 * f1^6 * f′0^5 * f′1^9 * f′′′0 -
        1680 * f0^8 * f1^6 * f′0^4 * f′1^9 * f′′0^2 +
        9900 * f0^8 * f1^5 * f′0^9 * f′1^6 * f′′1 +
        19800 * f0^8 * f1^4 * f′0^9 * f′1^8 - 19800 * f0^8 * f1^3 * f′0^9 * f′1^9 * x1 +
        5 * f0^7 * f1^9 * f′0^9 * f′1^3 * f′′′′′1 -
        75 * f0^7 * f1^9 * f′0^9 * f′1^2 * f′′1 * f′′′′1 -
        50 * f0^7 * f1^9 * f′0^9 * f′1^2 * f′′′1^2 +
        525 * f0^7 * f1^9 * f′0^9 * f′1 * f′′1^2 * f′′′1 -
        525 * f0^7 * f1^9 * f′0^9 * f′′1^4 - 10 * f0^7 * f1^9 * f′0^3 * f′1^9 * f′′′′′0 +
        150 * f0^7 * f1^9 * f′0^2 * f′1^9 * f′′0 * f′′′′0 +
        100 * f0^7 * f1^9 * f′0^2 * f′1^9 * f′′′0^2 -
        1050 * f0^7 * f1^9 * f′0 * f′1^9 * f′′0^2 * f′′′0 +
        1050 * f0^7 * f1^9 * f′1^9 * f′′0^4 - 145 * f0^7 * f1^8 * f′0^9 * f′1^4 * f′′′′1 -
        1450 * f0^7 * f1^8 * f′0^9 * f′1^3 * f′′1 * f′′′1 -
        2175 * f0^7 * f1^8 * f′0^9 * f′1^2 * f′′1^2 -
        230 * f0^7 * f1^8 * f′0^4 * f′1^9 * f′′′′0 -
        2300 * f0^7 * f1^8 * f′0^3 * f′1^9 * f′′0 * f′′′0 -
        3450 * f0^7 * f1^8 * f′0^2 * f′1^9 * f′′0^2 +
        1840 * f0^7 * f1^7 * f′0^9 * f′1^5 * f′′′1 -
        5520 * f0^7 * f1^7 * f′0^9 * f′1^4 * f′′1^2 -
        1840 * f0^7 * f1^7 * f′0^5 * f′1^9 * f′′′0 +
        5520 * f0^7 * f1^7 * f′0^4 * f′1^9 * f′′0^2 -
        12240 * f0^7 * f1^6 * f′0^9 * f′1^6 * f′′1 -
        5040 * f0^7 * f1^6 * f′0^6 * f′1^9 * f′′0 - 39600 * f0^7 * f1^5 * f′0^9 * f′1^8 +
        39600 * f0^7 * f1^4 * f′0^9 * f′1^9 * x1 - f0^6 * f1^10 * f′0^9 * f′1^3 * f′′′′′1 +
        15 * f0^6 * f1^10 * f′0^9 * f′1^2 * f′′1 * f′′′′1 +
        10 * f0^6 * f1^10 * f′0^9 * f′1^2 * f′′′1^2 -
        105 * f0^6 * f1^10 * f′0^9 * f′1 * f′′1^2 * f′′′1 +
        105 * f0^6 * f1^10 * f′0^9 * f′′1^4 +
        5 * f0^6 * f1^10 * f′0^3 * f′1^9 * f′′′′′0 -
        75 * f0^6 * f1^10 * f′0^2 * f′1^9 * f′′0 * f′′′′0 -
        50 * f0^6 * f1^10 * f′0^2 * f′1^9 * f′′′0^2 +
        525 * f0^6 * f1^10 * f′0 * f′1^9 * f′′0^2 * f′′′0 -
        525 * f0^6 * f1^10 * f′1^9 * f′′0^4 +
        35 * f0^6 * f1^9 * f′0^9 * f′1^4 * f′′′′1 +
        350 * f0^6 * f1^9 * f′0^9 * f′1^3 * f′′1 * f′′′1 +
        525 * f0^6 * f1^9 * f′0^9 * f′1^2 * f′′1^2 +
        170 * f0^6 * f1^9 * f′0^4 * f′1^9 * f′′′′0 +
        1700 * f0^6 * f1^9 * f′0^3 * f′1^9 * f′′0 * f′′′0 +
        2550 * f0^6 * f1^9 * f′0^2 * f′1^9 * f′′0^2 -
        560 * f0^6 * f1^8 * f′0^9 * f′1^5 * f′′′1 +
        1680 * f0^6 * f1^8 * f′0^9 * f′1^4 * f′′1^2 +
        2180 * f0^6 * f1^8 * f′0^5 * f′1^9 * f′′′0 -
        6540 * f0^6 * f1^8 * f′0^4 * f′1^9 * f′′0^2 +
        5040 * f0^6 * f1^7 * f′0^9 * f′1^6 * f′′1 +
        12240 * f0^6 * f1^7 * f′0^6 * f′1^9 * f′′0 +
        25200 * f0^6 * f1^6 * f′0^9 * f′1^8 - 25200 * f0^6 * f1^6 * f′0^8 * f′1^9 -
        55440 * f0^6 * f1^5 * f′0^9 * f′1^9 * x1 - f0^5 * f1^11 * f′0^3 * f′1^9 * f′′′′′0 +
        15 * f0^5 * f1^11 * f′0^2 * f′1^9 * f′′0 * f′′′′0 +
        10 * f0^5 * f1^11 * f′0^2 * f′1^9 * f′′′0^2 -
        105 * f0^5 * f1^11 * f′0 * f′1^9 * f′′0^2 * f′′′0 +
        105 * f0^5 * f1^11 * f′1^9 * f′′0^4 - 55 * f0^5 * f1^10 * f′0^4 * f′1^9 * f′′′′0 -
        550 * f0^5 * f1^10 * f′0^3 * f′1^9 * f′′0 * f′′′0 -
        825 * f0^5 * f1^10 * f′0^2 * f′1^9 * f′′0^2 -
        1100 * f0^5 * f1^9 * f′0^5 * f′1^9 * f′′′0 +
        3300 * f0^5 * f1^9 * f′0^4 * f′1^9 * f′′0^2 -
        9900 * f0^5 * f1^8 * f′0^6 * f′1^9 * f′′0 +
        39600 * f0^5 * f1^7 * f′0^8 * f′1^9 +
        55440 * f0^5 * f1^6 * f′0^9 * f′1^9 * x0 +
        5 * f0^4 * f1^11 * f′0^4 * f′1^9 * f′′′′0 +
        50 * f0^4 * f1^11 * f′0^3 * f′1^9 * f′′0 * f′′′0 +
        75 * f0^4 * f1^11 * f′0^2 * f′1^9 * f′′0^2 +
        220 * f0^4 * f1^10 * f′0^5 * f′1^9 * f′′′0 -
        660 * f0^4 * f1^10 * f′0^4 * f′1^9 * f′′0^2 +
        3300 * f0^4 * f1^9 * f′0^6 * f′1^9 * f′′0 - 19800 * f0^4 * f1^8 * f′0^8 * f′1^9 -
        39600 * f0^4 * f1^7 * f′0^9 * f′1^9 * x0 -
        20 * f0^3 * f1^11 * f′0^5 * f′1^9 * f′′′0 +
        60 * f0^3 * f1^11 * f′0^4 * f′1^9 * f′′0^2 -
        660 * f0^3 * f1^10 * f′0^6 * f′1^9 * f′′0 +
        6600 * f0^3 * f1^9 * f′0^8 * f′1^9 +
        19800 * f0^3 * f1^8 * f′0^9 * f′1^9 * x0 +
        60 * f0^2 * f1^11 * f′0^6 * f′1^9 * f′′0 - 1320 * f0^2 * f1^10 * f′0^8 * f′1^9 -
        6600 * f0^2 * f1^9 * f′0^9 * f′1^9 * x0 +
        120 * f0 * f1^11 * f′0^8 * f′1^9 +
        1320 * f0 * f1^10 * f′0^9 * f′1^9 * x0 - 120 * f1^11 * f′0^9 * f′1^9 * x0
    ) / (
        120 *
        f′0^9 *
        f′1^9 *
        (
            f0^11 - 11 * f0^10 * f1 + 55 * f0^9 * f1^2 - 165 * f0^8 * f1^3 +
            330 * f0^7 * f1^4 - 462 * f0^6 * f1^5 + 462 * f0^5 * f1^6 - 330 * f0^4 * f1^7 +
            165 * f0^3 * f1^8 - 55 * f0^2 * f1^9 + 11 * f0 * f1^10 - f1^11
        )
    )
end

function lmm(::LithBoonkkampIJzerman{3,5}, xs, fs, f′s, f′′s, f′′′s, f′′′′s, f′′′′′s)
    x0, x1, x2 = xs
    f0, f1, f2 = fs
    f′0, f′1, f′2 = f′s
    f′′0, f′′1, f′′2 = f′′s
    f′′′0, f′′′1, f′′′2 = f′′′s
    f′′′′0, f′′′′1, f′′′′2 = f′′′′s
    f′′′′′0, f′′′′′1, f′′′′′2 = f′′′′′s

    error("not computed")
end

## n = 6
function lmm(
    ::LithBoonkkampIJzerman{1,6},
    xs,
    fs,
    f′s,
    f′′s,
    f′′′s,
    f′′′′s,
    f′′′′′s,
    f′′′′′′s,
)
    x0 = xs[1]
    f0 = fs[1]
    f′0 = f′s[1]
    f′′0 = f′′s[1]
    f′′′0 = f′′′s[1]
    f′′′′0 = f′′′′s[1]
    f′′′′′0 = f′′′′′s[1]
    f′′′′′′0 = f′′′′′′s[1]

    (
        f0^6 * (
            -f′0^4 * f′′′′′′0 + 21 * f′0^3 * f′′0 * f′′′′′0 + 35 * f′0^3 * f′′′0 * f′′′′0 -
            210 * f′0^2 * f′′0^2 * f′′′′0 - 280 * f′0^2 * f′′0 * f′′′0^2 +
            1260 * f′0 * f′′0^3 * f′′′0 - 945 * f′′0^5
        ) +
        6 *
        f0^5 *
        f′0^2 *
        (
            f′0^3 * f′′′′′0 - 15 * f′0^2 * f′′0 * f′′′′0 - 10 * f′0^2 * f′′′0^2 +
            105 * f′0 * f′′0^2 * f′′′0 - 105 * f′′0^4
        ) - 30 * f0^4 * f′0^4 * (f′0^2 * f′′′′0 + 10 * f′0 * f′′0 * f′′′0 + 15 * f′′0^2) +
        120 * f0^3 * f′0^6 * (f′0 * f′′′0 - 3 * f′′0^2) - 360 * f0^2 * f′0^8 * f′′0 -
        720 * f0 * f′0^10 + 720 * f′0^11 * x0
    ) / (720 * f′0^11)
end

function lmm(
    ::LithBoonkkampIJzerman{2,6},
    xs,
    fs,
    f′s,
    f′′s,
    f′′′s,
    f′′′′s,
    f′′′′′s,
    f′′′′′′s,
)
    x0, x1 = xs
    f0, f1 = fs
    f′0, f′1 = f′s
    f′′0, f′′1 = f′′s
    f′′′0, f′′′1 = f′′′s
    f′′′′0, f′′′′1 = f′′′′s
    f′′′′′0, f′′′′′1 = f′′′′′s
    f′′′′′′0, f′′′′′′1 = f′′′′′′s

    (
        -f0^13 * f1^6 * f′0^11 * f′1^4 * f′′′′′′1 +
        21 * f0^13 * f1^6 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 +
        35 * f0^13 * f1^6 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 -
        210 * f0^13 * f1^6 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 -
        280 * f0^13 * f1^6 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 +
        1260 * f0^13 * f1^6 * f′0^11 * f′1 * f′′1^3 * f′′′1 -
        945 * f0^13 * f1^6 * f′0^11 * f′′1^5 + 6 * f0^13 * f1^5 * f′0^11 * f′1^5 * f′′′′′1 -
        90 * f0^13 * f1^5 * f′0^11 * f′1^4 * f′′1 * f′′′′1 -
        60 * f0^13 * f1^5 * f′0^11 * f′1^4 * f′′′1^2 +
        630 * f0^13 * f1^5 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 -
        630 * f0^13 * f1^5 * f′0^11 * f′1^2 * f′′1^4 -
        30 * f0^13 * f1^4 * f′0^11 * f′1^6 * f′′′′1 -
        300 * f0^13 * f1^4 * f′0^11 * f′1^5 * f′′1 * f′′′1 -
        450 * f0^13 * f1^4 * f′0^11 * f′1^4 * f′′1^2 +
        120 * f0^13 * f1^3 * f′0^11 * f′1^7 * f′′′1 -
        360 * f0^13 * f1^3 * f′0^11 * f′1^6 * f′′1^2 -
        360 * f0^13 * f1^2 * f′0^11 * f′1^8 * f′′1 - 720 * f0^13 * f1 * f′0^11 * f′1^10 +
        720 * f0^13 * f′0^11 * f′1^11 * x1 +
        6 * f0^12 * f1^7 * f′0^11 * f′1^4 * f′′′′′′1 -
        126 * f0^12 * f1^7 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 -
        210 * f0^12 * f1^7 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 +
        1260 * f0^12 * f1^7 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 +
        1680 * f0^12 * f1^7 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 -
        7560 * f0^12 * f1^7 * f′0^11 * f′1 * f′′1^3 * f′′′1 +
        5670 * f0^12 * f1^7 * f′0^11 * f′′1^5 +
        f0^12 * f1^7 * f′0^4 * f′1^11 * f′′′′′′0 -
        21 * f0^12 * f1^7 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 -
        35 * f0^12 * f1^7 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 +
        210 * f0^12 * f1^7 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 +
        280 * f0^12 * f1^7 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 -
        1260 * f0^12 * f1^7 * f′0 * f′1^11 * f′′0^3 * f′′′0 +
        945 * f0^12 * f1^7 * f′1^11 * f′′0^5 -
        78 * f0^12 * f1^6 * f′0^11 * f′1^5 * f′′′′′1 +
        1170 * f0^12 * f1^6 * f′0^11 * f′1^4 * f′′1 * f′′′′1 +
        780 * f0^12 * f1^6 * f′0^11 * f′1^4 * f′′′1^2 -
        8190 * f0^12 * f1^6 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 +
        8190 * f0^12 * f1^6 * f′0^11 * f′1^2 * f′′1^4 +
        390 * f0^12 * f1^5 * f′0^11 * f′1^6 * f′′′′1 +
        3900 * f0^12 * f1^5 * f′0^11 * f′1^5 * f′′1 * f′′′1 +
        5850 * f0^12 * f1^5 * f′0^11 * f′1^4 * f′′1^2 -
        1560 * f0^12 * f1^4 * f′0^11 * f′1^7 * f′′′1 +
        4680 * f0^12 * f1^4 * f′0^11 * f′1^6 * f′′1^2 +
        4680 * f0^12 * f1^3 * f′0^11 * f′1^8 * f′′1 +
        9360 * f0^12 * f1^2 * f′0^11 * f′1^10 - 9360 * f0^12 * f1 * f′0^11 * f′1^11 * x1 -
        15 * f0^11 * f1^8 * f′0^11 * f′1^4 * f′′′′′′1 +
        315 * f0^11 * f1^8 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 +
        525 * f0^11 * f1^8 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 -
        3150 * f0^11 * f1^8 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 -
        4200 * f0^11 * f1^8 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 +
        18900 * f0^11 * f1^8 * f′0^11 * f′1 * f′′1^3 * f′′′1 -
        14175 * f0^11 * f1^8 * f′0^11 * f′′1^5 -
        6 * f0^11 * f1^8 * f′0^4 * f′1^11 * f′′′′′′0 +
        126 * f0^11 * f1^8 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 +
        210 * f0^11 * f1^8 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 -
        1260 * f0^11 * f1^8 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 -
        1680 * f0^11 * f1^8 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 +
        7560 * f0^11 * f1^8 * f′0 * f′1^11 * f′′0^3 * f′′′0 -
        5670 * f0^11 * f1^8 * f′1^11 * f′′0^5 +
        300 * f0^11 * f1^7 * f′0^11 * f′1^5 * f′′′′′1 -
        4500 * f0^11 * f1^7 * f′0^11 * f′1^4 * f′′1 * f′′′′1 -
        3000 * f0^11 * f1^7 * f′0^11 * f′1^4 * f′′′1^2 +
        31500 * f0^11 * f1^7 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 -
        31500 * f0^11 * f1^7 * f′0^11 * f′1^2 * f′′1^4 -
        48 * f0^11 * f1^7 * f′0^5 * f′1^11 * f′′′′′0 +
        720 * f0^11 * f1^7 * f′0^4 * f′1^11 * f′′0 * f′′′′0 +
        480 * f0^11 * f1^7 * f′0^4 * f′1^11 * f′′′0^2 -
        5040 * f0^11 * f1^7 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 +
        5040 * f0^11 * f1^7 * f′0^2 * f′1^11 * f′′0^4 -
        2340 * f0^11 * f1^6 * f′0^11 * f′1^6 * f′′′′1 -
        23400 * f0^11 * f1^6 * f′0^11 * f′1^5 * f′′1 * f′′′1 -
        35100 * f0^11 * f1^6 * f′0^11 * f′1^4 * f′′1^2 +
        9360 * f0^11 * f1^5 * f′0^11 * f′1^7 * f′′′1 -
        28080 * f0^11 * f1^5 * f′0^11 * f′1^6 * f′′1^2 -
        28080 * f0^11 * f1^4 * f′0^11 * f′1^8 * f′′1 -
        56160 * f0^11 * f1^3 * f′0^11 * f′1^10 +
        56160 * f0^11 * f1^2 * f′0^11 * f′1^11 * x1 +
        20 * f0^10 * f1^9 * f′0^11 * f′1^4 * f′′′′′′1 -
        420 * f0^10 * f1^9 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 -
        700 * f0^10 * f1^9 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 +
        4200 * f0^10 * f1^9 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 +
        5600 * f0^10 * f1^9 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 -
        25200 * f0^10 * f1^9 * f′0^11 * f′1 * f′′1^3 * f′′′1 +
        18900 * f0^10 * f1^9 * f′0^11 * f′′1^5 +
        15 * f0^10 * f1^9 * f′0^4 * f′1^11 * f′′′′′′0 -
        315 * f0^10 * f1^9 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 -
        525 * f0^10 * f1^9 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 +
        3150 * f0^10 * f1^9 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 +
        4200 * f0^10 * f1^9 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 -
        18900 * f0^10 * f1^9 * f′0 * f′1^11 * f′′0^3 * f′′′0 +
        14175 * f0^10 * f1^9 * f′1^11 * f′′0^5 -
        540 * f0^10 * f1^8 * f′0^11 * f′1^5 * f′′′′′1 +
        8100 * f0^10 * f1^8 * f′0^11 * f′1^4 * f′′1 * f′′′′1 +
        5400 * f0^10 * f1^8 * f′0^11 * f′1^4 * f′′′1^2 -
        56700 * f0^10 * f1^8 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 +
        56700 * f0^10 * f1^8 * f′0^11 * f′1^2 * f′′1^4 +
        246 * f0^10 * f1^8 * f′0^5 * f′1^11 * f′′′′′0 -
        3690 * f0^10 * f1^8 * f′0^4 * f′1^11 * f′′0 * f′′′′0 -
        2460 * f0^10 * f1^8 * f′0^4 * f′1^11 * f′′′0^2 +
        25830 * f0^10 * f1^8 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 -
        25830 * f0^10 * f1^8 * f′0^2 * f′1^11 * f′′0^4 +
        6060 * f0^10 * f1^7 * f′0^11 * f′1^6 * f′′′′1 +
        60600 * f0^10 * f1^7 * f′0^11 * f′1^5 * f′′1 * f′′′1 +
        90900 * f0^10 * f1^7 * f′0^11 * f′1^4 * f′′1^2 +
        1080 * f0^10 * f1^7 * f′0^6 * f′1^11 * f′′′′0 +
        10800 * f0^10 * f1^7 * f′0^5 * f′1^11 * f′′0 * f′′′0 +
        16200 * f0^10 * f1^7 * f′0^4 * f′1^11 * f′′0^2 -
        34320 * f0^10 * f1^6 * f′0^11 * f′1^7 * f′′′1 +
        102960 * f0^10 * f1^6 * f′0^11 * f′1^6 * f′′1^2 +
        102960 * f0^10 * f1^5 * f′0^11 * f′1^8 * f′′1 +
        205920 * f0^10 * f1^4 * f′0^11 * f′1^10 -
        205920 * f0^10 * f1^3 * f′0^11 * f′1^11 * x1 -
        15 * f0^9 * f1^10 * f′0^11 * f′1^4 * f′′′′′′1 +
        315 * f0^9 * f1^10 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 +
        525 * f0^9 * f1^10 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 -
        3150 * f0^9 * f1^10 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 -
        4200 * f0^9 * f1^10 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 +
        18900 * f0^9 * f1^10 * f′0^11 * f′1 * f′′1^3 * f′′′1 -
        14175 * f0^9 * f1^10 * f′0^11 * f′′1^5 -
        20 * f0^9 * f1^10 * f′0^4 * f′1^11 * f′′′′′′0 +
        420 * f0^9 * f1^10 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 +
        700 * f0^9 * f1^10 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 -
        4200 * f0^9 * f1^10 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 -
        5600 * f0^9 * f1^10 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 +
        25200 * f0^9 * f1^10 * f′0 * f′1^11 * f′′0^3 * f′′′0 -
        18900 * f0^9 * f1^10 * f′1^11 * f′′0^5 +
        510 * f0^9 * f1^9 * f′0^11 * f′1^5 * f′′′′′1 -
        7650 * f0^9 * f1^9 * f′0^11 * f′1^4 * f′′1 * f′′′′1 -
        5100 * f0^9 * f1^9 * f′0^11 * f′1^4 * f′′′1^2 +
        53550 * f0^9 * f1^9 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 -
        53550 * f0^9 * f1^9 * f′0^11 * f′1^2 * f′′1^4 -
        510 * f0^9 * f1^9 * f′0^5 * f′1^11 * f′′′′′0 +
        7650 * f0^9 * f1^9 * f′0^4 * f′1^11 * f′′0 * f′′′′0 +
        5100 * f0^9 * f1^9 * f′0^4 * f′1^11 * f′′′0^2 -
        53550 * f0^9 * f1^9 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 +
        53550 * f0^9 * f1^9 * f′0^2 * f′1^11 * f′′0^4 -
        7590 * f0^9 * f1^8 * f′0^11 * f′1^6 * f′′′′1 -
        75900 * f0^9 * f1^8 * f′0^11 * f′1^5 * f′′1 * f′′′1 -
        113850 * f0^9 * f1^8 * f′0^11 * f′1^4 * f′′1^2 -
        4590 * f0^9 * f1^8 * f′0^6 * f′1^11 * f′′′′0 -
        45900 * f0^9 * f1^8 * f′0^5 * f′1^11 * f′′0 * f′′′0 -
        68850 * f0^9 * f1^8 * f′0^4 * f′1^11 * f′′0^2 +
        60600 * f0^9 * f1^7 * f′0^11 * f′1^7 * f′′′1 -
        181800 * f0^9 * f1^7 * f′0^11 * f′1^6 * f′′1^2 -
        14400 * f0^9 * f1^7 * f′0^7 * f′1^11 * f′′′0 +
        43200 * f0^9 * f1^7 * f′0^6 * f′1^11 * f′′0^2 -
        257400 * f0^9 * f1^6 * f′0^11 * f′1^8 * f′′1 -
        514800 * f0^9 * f1^5 * f′0^11 * f′1^10 +
        514800 * f0^9 * f1^4 * f′0^11 * f′1^11 * x1 +
        6 * f0^8 * f1^11 * f′0^11 * f′1^4 * f′′′′′′1 -
        126 * f0^8 * f1^11 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 -
        210 * f0^8 * f1^11 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 +
        1260 * f0^8 * f1^11 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 +
        1680 * f0^8 * f1^11 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 -
        7560 * f0^8 * f1^11 * f′0^11 * f′1 * f′′1^3 * f′′′1 +
        5670 * f0^8 * f1^11 * f′0^11 * f′′1^5 +
        15 * f0^8 * f1^11 * f′0^4 * f′1^11 * f′′′′′′0 -
        315 * f0^8 * f1^11 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 -
        525 * f0^8 * f1^11 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 +
        3150 * f0^8 * f1^11 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 +
        4200 * f0^8 * f1^11 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 -
        18900 * f0^8 * f1^11 * f′0 * f′1^11 * f′′0^3 * f′′′0 +
        14175 * f0^8 * f1^11 * f′1^11 * f′′0^5 -
        246 * f0^8 * f1^10 * f′0^11 * f′1^5 * f′′′′′1 +
        3690 * f0^8 * f1^10 * f′0^11 * f′1^4 * f′′1 * f′′′′1 +
        2460 * f0^8 * f1^10 * f′0^11 * f′1^4 * f′′′1^2 -
        25830 * f0^8 * f1^10 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 +
        25830 * f0^8 * f1^10 * f′0^11 * f′1^2 * f′′1^4 +
        540 * f0^8 * f1^10 * f′0^5 * f′1^11 * f′′′′′0 -
        8100 * f0^8 * f1^10 * f′0^4 * f′1^11 * f′′0 * f′′′′0 -
        5400 * f0^8 * f1^10 * f′0^4 * f′1^11 * f′′′0^2 +
        56700 * f0^8 * f1^10 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 -
        56700 * f0^8 * f1^10 * f′0^2 * f′1^11 * f′′0^4 +
        4590 * f0^8 * f1^9 * f′0^11 * f′1^6 * f′′′′1 +
        45900 * f0^8 * f1^9 * f′0^11 * f′1^5 * f′′1 * f′′′1 +
        68850 * f0^8 * f1^9 * f′0^11 * f′1^4 * f′′1^2 +
        7590 * f0^8 * f1^9 * f′0^6 * f′1^11 * f′′′′0 +
        75900 * f0^8 * f1^9 * f′0^5 * f′1^11 * f′′0 * f′′′0 +
        113850 * f0^8 * f1^9 * f′0^4 * f′1^11 * f′′0^2 -
        48600 * f0^8 * f1^8 * f′0^11 * f′1^7 * f′′′1 +
        145800 * f0^8 * f1^8 * f′0^11 * f′1^6 * f′′1^2 +
        48600 * f0^8 * f1^8 * f′0^7 * f′1^11 * f′′′0 -
        145800 * f0^8 * f1^8 * f′0^6 * f′1^11 * f′′0^2 +
        297000 * f0^8 * f1^7 * f′0^11 * f′1^8 * f′′1 +
        118800 * f0^8 * f1^7 * f′0^8 * f′1^11 * f′′0 +
        926640 * f0^8 * f1^6 * f′0^11 * f′1^10 -
        926640 * f0^8 * f1^5 * f′0^11 * f′1^11 * x1 -
        f0^7 * f1^12 * f′0^11 * f′1^4 * f′′′′′′1 +
        21 * f0^7 * f1^12 * f′0^11 * f′1^3 * f′′1 * f′′′′′1 +
        35 * f0^7 * f1^12 * f′0^11 * f′1^3 * f′′′1 * f′′′′1 -
        210 * f0^7 * f1^12 * f′0^11 * f′1^2 * f′′1^2 * f′′′′1 -
        280 * f0^7 * f1^12 * f′0^11 * f′1^2 * f′′1 * f′′′1^2 +
        1260 * f0^7 * f1^12 * f′0^11 * f′1 * f′′1^3 * f′′′1 -
        945 * f0^7 * f1^12 * f′0^11 * f′′1^5 -
        6 * f0^7 * f1^12 * f′0^4 * f′1^11 * f′′′′′′0 +
        126 * f0^7 * f1^12 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 +
        210 * f0^7 * f1^12 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 -
        1260 * f0^7 * f1^12 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 -
        1680 * f0^7 * f1^12 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 +
        7560 * f0^7 * f1^12 * f′0 * f′1^11 * f′′0^3 * f′′′0 -
        5670 * f0^7 * f1^12 * f′1^11 * f′′0^5 +
        48 * f0^7 * f1^11 * f′0^11 * f′1^5 * f′′′′′1 -
        720 * f0^7 * f1^11 * f′0^11 * f′1^4 * f′′1 * f′′′′1 -
        480 * f0^7 * f1^11 * f′0^11 * f′1^4 * f′′′1^2 +
        5040 * f0^7 * f1^11 * f′0^11 * f′1^3 * f′′1^2 * f′′′1 -
        5040 * f0^7 * f1^11 * f′0^11 * f′1^2 * f′′1^4 -
        300 * f0^7 * f1^11 * f′0^5 * f′1^11 * f′′′′′0 +
        4500 * f0^7 * f1^11 * f′0^4 * f′1^11 * f′′0 * f′′′′0 +
        3000 * f0^7 * f1^11 * f′0^4 * f′1^11 * f′′′0^2 -
        31500 * f0^7 * f1^11 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 +
        31500 * f0^7 * f1^11 * f′0^2 * f′1^11 * f′′0^4 -
        1080 * f0^7 * f1^10 * f′0^11 * f′1^6 * f′′′′1 -
        10800 * f0^7 * f1^10 * f′0^11 * f′1^5 * f′′1 * f′′′1 -
        16200 * f0^7 * f1^10 * f′0^11 * f′1^4 * f′′1^2 -
        6060 * f0^7 * f1^10 * f′0^6 * f′1^11 * f′′′′0 -
        60600 * f0^7 * f1^10 * f′0^5 * f′1^11 * f′′0 * f′′′0 -
        90900 * f0^7 * f1^10 * f′0^4 * f′1^11 * f′′0^2 +
        14400 * f0^7 * f1^9 * f′0^11 * f′1^7 * f′′′1 -
        43200 * f0^7 * f1^9 * f′0^11 * f′1^6 * f′′1^2 -
        60600 * f0^7 * f1^9 * f′0^7 * f′1^11 * f′′′0 +
        181800 * f0^7 * f1^9 * f′0^6 * f′1^11 * f′′0^2 -
        118800 * f0^7 * f1^8 * f′0^11 * f′1^8 * f′′1 -
        297000 * f0^7 * f1^8 * f′0^8 * f′1^11 * f′′0 -
        570240 * f0^7 * f1^7 * f′0^11 * f′1^10 +
        570240 * f0^7 * f1^7 * f′0^10 * f′1^11 +
        1235520 * f0^7 * f1^6 * f′0^11 * f′1^11 * x1 +
        f0^6 * f1^13 * f′0^4 * f′1^11 * f′′′′′′0 -
        21 * f0^6 * f1^13 * f′0^3 * f′1^11 * f′′0 * f′′′′′0 -
        35 * f0^6 * f1^13 * f′0^3 * f′1^11 * f′′′0 * f′′′′0 +
        210 * f0^6 * f1^13 * f′0^2 * f′1^11 * f′′0^2 * f′′′′0 +
        280 * f0^6 * f1^13 * f′0^2 * f′1^11 * f′′0 * f′′′0^2 -
        1260 * f0^6 * f1^13 * f′0 * f′1^11 * f′′0^3 * f′′′0 +
        945 * f0^6 * f1^13 * f′1^11 * f′′0^5 +
        78 * f0^6 * f1^12 * f′0^5 * f′1^11 * f′′′′′0 -
        1170 * f0^6 * f1^12 * f′0^4 * f′1^11 * f′′0 * f′′′′0 -
        780 * f0^6 * f1^12 * f′0^4 * f′1^11 * f′′′0^2 +
        8190 * f0^6 * f1^12 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 -
        8190 * f0^6 * f1^12 * f′0^2 * f′1^11 * f′′0^4 +
        2340 * f0^6 * f1^11 * f′0^6 * f′1^11 * f′′′′0 +
        23400 * f0^6 * f1^11 * f′0^5 * f′1^11 * f′′0 * f′′′0 +
        35100 * f0^6 * f1^11 * f′0^4 * f′1^11 * f′′0^2 +
        34320 * f0^6 * f1^10 * f′0^7 * f′1^11 * f′′′0 -
        102960 * f0^6 * f1^10 * f′0^6 * f′1^11 * f′′0^2 +
        257400 * f0^6 * f1^9 * f′0^8 * f′1^11 * f′′0 -
        926640 * f0^6 * f1^8 * f′0^10 * f′1^11 -
        1235520 * f0^6 * f1^7 * f′0^11 * f′1^11 * x0 -
        6 * f0^5 * f1^13 * f′0^5 * f′1^11 * f′′′′′0 +
        90 * f0^5 * f1^13 * f′0^4 * f′1^11 * f′′0 * f′′′′0 +
        60 * f0^5 * f1^13 * f′0^4 * f′1^11 * f′′′0^2 -
        630 * f0^5 * f1^13 * f′0^3 * f′1^11 * f′′0^2 * f′′′0 +
        630 * f0^5 * f1^13 * f′0^2 * f′1^11 * f′′0^4 -
        390 * f0^5 * f1^12 * f′0^6 * f′1^11 * f′′′′0 -
        3900 * f0^5 * f1^12 * f′0^5 * f′1^11 * f′′0 * f′′′0 -
        5850 * f0^5 * f1^12 * f′0^4 * f′1^11 * f′′0^2 -
        9360 * f0^5 * f1^11 * f′0^7 * f′1^11 * f′′′0 +
        28080 * f0^5 * f1^11 * f′0^6 * f′1^11 * f′′0^2 -
        102960 * f0^5 * f1^10 * f′0^8 * f′1^11 * f′′0 +
        514800 * f0^5 * f1^9 * f′0^10 * f′1^11 +
        926640 * f0^5 * f1^8 * f′0^11 * f′1^11 * x0 +
        30 * f0^4 * f1^13 * f′0^6 * f′1^11 * f′′′′0 +
        300 * f0^4 * f1^13 * f′0^5 * f′1^11 * f′′0 * f′′′0 +
        450 * f0^4 * f1^13 * f′0^4 * f′1^11 * f′′0^2 +
        1560 * f0^4 * f1^12 * f′0^7 * f′1^11 * f′′′0 -
        4680 * f0^4 * f1^12 * f′0^6 * f′1^11 * f′′0^2 +
        28080 * f0^4 * f1^11 * f′0^8 * f′1^11 * f′′0 -
        205920 * f0^4 * f1^10 * f′0^10 * f′1^11 -
        514800 * f0^4 * f1^9 * f′0^11 * f′1^11 * x0 -
        120 * f0^3 * f1^13 * f′0^7 * f′1^11 * f′′′0 +
        360 * f0^3 * f1^13 * f′0^6 * f′1^11 * f′′0^2 -
        4680 * f0^3 * f1^12 * f′0^8 * f′1^11 * f′′0 +
        56160 * f0^3 * f1^11 * f′0^10 * f′1^11 +
        205920 * f0^3 * f1^10 * f′0^11 * f′1^11 * x0 +
        360 * f0^2 * f1^13 * f′0^8 * f′1^11 * f′′0 - 9360 * f0^2 * f1^12 * f′0^10 * f′1^11 -
        56160 * f0^2 * f1^11 * f′0^11 * f′1^11 * x0 +
        720 * f0 * f1^13 * f′0^10 * f′1^11 +
        9360 * f0 * f1^12 * f′0^11 * f′1^11 * x0 - 720 * f1^13 * f′0^11 * f′1^11 * x0
    ) / (
        720 *
        f′0^11 *
        f′1^11 *
        (
            f0^13 - 13 * f0^12 * f1 + 78 * f0^11 * f1^2 - 286 * f0^10 * f1^3 +
            715 * f0^9 * f1^4 - 1287 * f0^8 * f1^5 + 1716 * f0^7 * f1^6 -
            1716 * f0^6 * f1^7 + 1287 * f0^5 * f1^8 - 715 * f0^4 * f1^9 +
            286 * f0^3 * f1^10 - 78 * f0^2 * f1^11 + 13 * f0 * f1^12 - f1^13
        )
    )
end

function lmm(
    ::LithBoonkkampIJzerman{3,6},
    xs,
    fs,
    f′s,
    f′′s,
    f′′′s,
    f′′′′s,
    f′′′′′s,
    f′′′′′′s,
)
    error("not computed")
end
