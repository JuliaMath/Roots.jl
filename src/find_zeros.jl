# Attempt to find all zeros in an interval (a,b)

# Algorithm due to @djsegal in https://github.com/JuliaMath/Roots.jl/pull/113

## -----
## Bisection
##
## Essentially from Jason Merrill https://gist.github.com/jwmerrill/9012954
## cf. http://squishythinking.com/2014/02/22/bisecting-floats/
## This also borrows a trick from https://discourse.julialang.org/t/simple-and-fast-bisection/14886
## where we keep x1 so that y1 is negative, and x2 so that y2 is positive
## this allows the use of signbit over y1*y2 < 0 which avoid < and a multiplication
## this has a small, but noticeable impact on performance.
"""
    bisection(f, a, b; [xatol, xrtol])

Performs bisection method to find a zero of a continuous
function.

It is assumed that `(a,b)` is a bracket, that is, the function has
different signs at `a` and `b`. The interval `(a,b)` is converted to floating point
and shrunk when `a` or `b` is infinite. The function `f` may be infinite for
the typical case. If `f` is not continuous, the algorithm may find
jumping points over the x axis, not just zeros.


If non-trivial tolerances are specified, the process will terminate
when the bracket `(a,b)` satisfies `isapprox(a, b, atol=xatol,
rtol=xrtol)`. For zero tolerances, the default, for `Float64`, `Float32`,
or `Float16` values, the process will terminate at a value `x` with
`f(x)=0` or `f(x)*f(prevfloat(x)) < 0 ` or `f(x) * f(nextfloat(x)) <
0`. For other number types, the `Roots.A42` method is used.

"""
function bisection(f, a::Number, b::Number; xatol=nothing, xrtol=nothing)
    x1, x2 = adjust_bracket(float.((a, b)))
    T = eltype(x1)

    atol = xatol === nothing ? zero(T) : abs(xatol)
    rtol = xrtol === nothing ? zero(one(T)) : abs(xrtol)
    CT = iszero(atol) && iszero(rtol) ? Val(:exact) : Val(:inexact)

    x1, x2 = float(x1), float(x2)
    y1, y2 = f(x1), f(x2)

    _unitless(y1 * y2) >= 0 && error("the interval provided does not bracket a root")

    if isneg(y2)
        x1, x2, y1, y2 = x2, x1, y2, y1
    end

    xm = Roots._middle(x1, x2) # for possibly mixed sign x1, x2
    ym = f(xm)

    while true
        if has_converged(CT, x1, x2, xm, ym, atol, rtol)
            return xm
        end

        if isneg(ym)
            x1, y1 = xm, ym
        else
            x2, y2 = xm, ym
        end

        xm = Roots.__middle(x1, x2)
        ym = f(xm)
    end
end

# -0.0 not returned by __middle, so isneg true on [-Inf, 0.0)
@inline isneg(x::T) where {T<:AbstractFloat} = signbit(x)
@inline isneg(x) = _unitless(x) < 0

@inline function has_converged(::Val{:exact}, x1, x2, m, ym, atol, rtol)
    iszero(ym) && return true
    isnan(ym) && return true
    x1 != m && m != x2 && return false
    return true
end

@inline function has_converged(::Val{:inexact}, x1, x2, m, ym, atol, rtol)
    iszero(ym) && return true
    isnan(ym) && return true
    val = abs(x1 - x2) <= atol + max(abs(x1), abs(x2)) * rtol
    return val
end

## ----
## --------------------------------------------------

## This is basically Order0(), but with different, default, tolerances employed
## It takes more function calls, but works harder to find exact zeros
## where exact means either iszero(fx), adjacent floats have sign change, or
## abs(fxn) <= 8 eps(xn)
"""
    dfree(f, xs)

A more robust secant method implementation

Solve for `f(x) = 0` using an algorithm from *Personal Calculator Has Key
to Solve Any Equation f(x) = 0*, the SOLVE button from the
[HP-34C](http://www.hpl.hp.com/hpjournal/pdfs/IssuePDFs/1979-12.pdf).

This is also implemented as the `Order0` method for `find_zero`.

The initial values can be specified as a pair of two values, as in
`(a,b)` or `[a,b]`, or as a single value, in which case a value of `b`
is computed, possibly from `fb`.  The basic idea is to follow the
secant method to convergence unless:

* a bracket is found, in which case `AlefeldPotraShi` is used;

* the secant method is not converging, in which case a few steps of a
  quadratic method are used to see if that improves matters.

Convergence occurs when `f(m) == 0`, there is a sign change between
`m` and an adjacent floating point value, or `f(m) <= 2^3*eps(m)`.

A value of `NaN` is returned if the algorithm takes too many steps
before identifying a zero.

# Examples

```julia
Roots.dfree(x -> x^5 - x - 1, 1.0)
```

"""
function dfree(f, xs)
    if length(xs) == 1
        a = float(xs[1])
        fa = f(a)

        h = eps(one(a))^(1 / 3)
        da = h * oneunit(a) + abs(a) * h^2 # adjust for if eps(a) > h
        b = float(a + da)
        fb = f(b)
    else
        a, b = promote(float(xs[1]), float(xs[2]))
        fa, fb = f(a), f(b)
    end

    nan = (0 * a) / (0 * a) # try to preserve type
    cnt, MAXCNT = 0, 5 * ceil(Int, -log(eps(one(a))))  # must be higher for BigFloat
    MAXQUAD = 3

    if abs(fa) > abs(fb)
        a, fa, b, fb = b, fb, a, fa
    end

    # we keep a, b, fa, fb, gamma, fgamma
    quad_ctr = 0
    while !iszero(fb)
        cnt += 1

        if sign(fa) * sign(fb) < 0
            return solve(ZeroProblem(f, (a, b))) # faster than bisection(f, a, b)
        end

        # take a secant step
        gamma = float(b - (b - a) * fb / (fb - fa))
        # modify if gamma is too small or too big
        if iszero(abs(gamma - b))
            gamma = b + 1 / 1000 * abs(b - a)  # too small
        elseif abs(gamma - b) >= 100 * abs(b - a)
            gamma = b + sign(gamma - b) * 100 * abs(b - a)  ## too big
        end
        fgamma = f(gamma)

        # change sign
        if sign(fgamma) * sign(fb) < 0
            return solve(ZeroProblem(f, (gamma, b))) # faster than bisection(f, gamma, b)
        end

        # decreasing
        if abs(fgamma) < abs(fb)
            a, fa, b, fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        gamma = float(quad_vertex(a, fa, b, fb, gamma, fgamma))
        fgamma = f(gamma)
        # decreasing now?
        if abs(fgamma) < abs(fb)
            a, fa, b, fb = b, fb, gamma, fgamma
            quad_ctr = 0
            cnt < MAXCNT && continue
        end

        quad_ctr += 1
        if (quad_ctr > MAXQUAD) || (cnt > MAXCNT) || iszero(gamma - b) || isnan(gamma)
            bprev, bnext = prevfloat(b), nextfloat(b)
            fbprev, fbnext = f(bprev), f(bnext)
            sign(fb) * sign(fbprev) < 0 && return b
            sign(fb) * sign(fbnext) < 0 && return b
            for (u, fu) in ((b, fb), (bprev, fbprev), (bnext, fbnext))
                abs(fu) / oneunit(fu) <= 2^3 * eps(u / oneunit(u)) && return u
            end
            return nan # Failed.
        end

        if abs(fgamma) < abs(fb)
            b, fb, a, fa = gamma, fgamma, b, fb
        else
            a, fa = gamma, fgamma
        end
    end
    b
end


# A naive approach to find zeros: split (a,b) by n points, look into each for a zero
# * k is oversampling rate for bisection. (It is relatively cheap to check for a bracket so we
#   oversample our intervals looking for brackets
# * assumes f(a) *not* a zero
function _fz(f, a, b, no_pts, k=4)
    zs = Real[]
    _fz!(zs, f, a, b, no_pts, k)
    zs
end

function _fz!(zs, f, a::T, b, no_pts, k::Int=4) where {T}
    pts = range(a, stop=b, length=(no_pts - 1) * k + 1)
    n::Int = length(pts)

    fs = f.(pts)
    sfs = sign.(fs)

    u::T = first(pts)  # keep track of bigger interval
    found_bisection_zero = false

    for (i, x) in enumerate(pts[1:end])
        q, r = divrem(i - 1, k)

        if i > 1 && iszero(r)
            v::T = x
            if !found_bisection_zero
                try
                    p1::T = identify_starting_point(u, v, sfs[(i - k):i])
                    rt::T = dfree(f, p1)
                    if !isnan(rt) && u < rt <= v
                        push!(zs, rt)
                    end
                catch err
                end
            end
            u = v
            found_bisection_zero = false
        end

        if i < n
            if iszero(fs[i + 1])
                found_bisection_zero = true # kinda
                push!(zs, pts[i + 1])
            elseif sfs[i] * sfs[i + 1] < 0
                found_bisection_zero = true
                rt = bisection(f, x, pts[i + 1])
                push!(zs, rt)
            end
        end
    end

    sort!(zs)
end

# the algorithm first scans for zeros using the naive approach, then
# splits (a,b) by these zeros. This struct holds the subintervals
struct Interval{T}
    a::T
    b::T
    depth::Int
end

Base.show(io::IO, alpha::Interval) = print(io, "($(alpha.a), $(alpha.b))")

# check if f(a) is non zero using tolerances max(atol, eps()), rtol
function _non_zero(fa, a::T, atol, rtol) where {T}
    a, r = atol, abs(a) * rtol * oneunit(fa) / oneunit(a), oneunit(fa) * eps(T)
    return abs(fa) >= max(promote(a, r)...)
end

# After splitting by zeros we have intervals (zm, zn) this is used to shrink
# to (zm+, zn-) where both are non-zeros, as defined above
function find_non_zero(f, a::T, barrier, xatol, xrtol, atol, rtol) where {T}
    nan = (0 * a) / (0 * a) # try to get typed NaN
    xtol = max(xatol, abs(a) * xrtol, oneunit(a) * eps(T))
    sgn = barrier > a ? 1 : -1
    ctr = 0
    x = a + 2^ctr * sgn * xtol
    while !_non_zero(f(x), x, atol, rtol)
        ctr += 1
        x += 2^ctr * sgn * xtol
        ((sgn > 0 && x > barrier) || (sgn < 0 && x < barrier)) && return nan
        ctr > 100 && return nan
    end

    x
end

# split a < z1 < z2 < ... < zn < b into intervals (a+,z1-), (z1+, z2-), ...
# where possible; push! onto ints
function make_intervals!(ints, f, a, b, zs, depth, xatol, xrtol, atol, rtol)
    pts = vcat(a, zs, b)

    for (u, v) in zip(pts[1:(end - 1)], pts[2:end])
        ur = find_non_zero(f, u, v, xatol, xrtol, atol, rtol)
        isnan(ur) && continue

        vl = find_non_zero(f, v, u, xatol, xrtol, atol, rtol)
        isnan(vl) && continue

        push!(ints, Interval(ur, vl, depth))
    end
end

# adjust what we mean by x1 ~ x2 for purposes of adding a new zero
function approx_close(z1, z2, xatol, xrtol)
    z₁, z₂, δ, ϵ = _unitless.((z1, z2, xatol, xrtol))
    return isapprox(z₁, z₂; atol=sqrt(δ), rtol=sqrt(ϵ))
end

# is proposed not near xs? (false and we add proposed)
function not_near(proposed, xs, xatol, xrtol)
    n = length(xs)
    n <= 1 && return true
    ind = n + 1
    for (i, rt) in enumerate(xs)
        if proposed < rt
            ind = i
            break
        end
    end
    if ind > 1 # xs[ind-1] <= propose < xs[ind]
        rt = xs[ind - 1]
        approx_close(proposed, rt, xatol, xrtol) && return false
    end
    if ind <= n # value to right
        rt = xs[ind]
        approx_close(proposed, rt, xatol, xrtol) && return false
    end
    return true
end

"""
    find_zeros(f, a, [b]; [no_pts=12, k=8, naive=false, xatol, xrtol, atol, rtol])

Search for zeros of `f` in the interval `[a,b]` with an heuristic
algorithm.


* `f`: a function or callable object
* `a`, `b`:  If `b` is specified, the interval ``[a,b]`` is used. If only `a` is specified, it is passed to `extrema` to define the interval to search over.   It is assumed that neither endpoint is a zero.

Returns a vector of zeros in sorted order, possibly empty.

# Extended help

# Examples

```jldoctest find_zeros
julia> using Roots

julia> find_zeros(x -> exp(x) - x^4, -5, 20)        # a few well-spaced zeros
3-element Vector{Float64}:
 -0.8155534188089606
  1.4296118247255556
  8.613169456441398

julia> find_zeros(x -> sin(x^2) + cos(x)^2, 0, 2pi)  # many zeros
12-element Vector{Float64}:
 1.78518032659534
 2.391345462376604
 3.2852368649448853
 3.3625557095737544
 4.016412952618305
 4.325091924521049
 4.68952781386834
 5.00494459113514
 5.35145266881871
 5.552319796014526
 5.974560835055425
 6.039177477770888

julia> find_zeros(x -> cos(x) + cos(2x), (0, 4pi))    # mix of simple, non-simple zeros
6-element Vector{Float64}:
  1.0471975511965976
  3.141592653589943
  5.235987755982988
  7.330382858376184
  9.424777960769228
 11.519173063162574

julia> f(x) = (x-0.5) * (x-0.5001) * (x-1)          # nearby zeros
f (generic function with 1 method)

julia> find_zeros(f, 0, 2)
3-element Vector{Float64}:
 0.5
 0.5001
 1.0

julia> f(x) = (x-0.5) * (x-0.5001) * (x-4) * (x-4.001) * (x-4.2)
f (generic function with 1 method)

julia> find_zeros(f, 0, 10)
3-element Vector{Float64}:
 0.5
 0.5001
 4.2

julia> f(x) = (x-0.5)^2 * (x-0.5001)^3 * (x-4) * (x-4.001) * (x-4.2)^2  # hard to identify
f (generic function with 1 method)

julia> find_zeros(f, 0, 10, no_pts=21)                # too hard for default
5-element Vector{Float64}:
 0.49999999999999994
 0.5001
 4.0
 4.001
 4.200000000000001
```

!!! note
    Some cases where the number of zeros may be underreported:
    * if the initial interval, `(a,b)`, is too wide
    * if there are zeros  that are very nearby
    * the function is flat, e.g., `x->0`.

----

The basic algorithm checks for zeros among the endpoints, and then
divides the interval `(a,b)` into `no_pts-1` subintervals and then
proceeds to look for zeros through bisection or a derivative-free
method.  As checking for a bracketing interval is relatively cheap and
bisection is guaranteed to converge, each interval has `k` pairs of
intermediate points checked for a bracket.

If any zeros are found, the algorithm uses these to partition `(a,b)`
into subintervals. Each subinterval is shrunk so that the endpoints
are not zeros and the process is repeated on the subinterval. If the
initial interval is too large, then the naive scanning for zeros may
be fruitless and no zeros will be reported. If there are nearby zeros,
the shrinking of the interval may jump over them, though as seen in
the examples, nearby roots can be identified correctly, though for
really nearby points, or very flat functions, it may help to increase
`no_pts`.


The tolerances are used to shrink the intervals, but not to find zeros
within a search. For searches, bisection is guaranteed to converge
with no specified tolerance. For the derivative free search, a
modification of the `Order0` method is used, which at worst case
compares `|f(x)| <= 8*eps(x)` to identify a zero. The algorithm might
identify more than one value for a zero, due to floating point
approximations. If a potential pair of zeros satisfy
`isapprox(a,b,atol=sqrt(xatol), rtol=sqrt(xrtol))` then they are
consolidated.

The algorithm can make many function calls. When zeros are found in an
interval, the naive search is carried out on each subinterval. To cut
down on function calls, though with some increased chance of missing
some zeros, the adaptive nature can be skipped with the argument
`naive=true` or the number of points stepped down.


The algorithm is derived from one in a
[PR](https://github.com/JuliaMath/Roots.jl/pull/113) by @djsegal.


!!! note
    The `IntervalRootFinding` package provides a rigorous alternative to this heuristic one.
    That package uses interval arithmetic, so can compute bounds on the size of the image of
    an interval under `f`. If this image includes `0`, then it can look for the zero.
    Bisection, on the other hand, only will look for a zero if the two endpoints have different signs,
    a much more rigid condition for a potential zero.

!!! note "`IntervalRootFinding` extension"
    As of version `1.9` an extension is provided so that when the `IntervalRootFinding` package is loaded,
    the `find_zeros` function will call `IntervalRootFinding.roots` to find the isolating brackets and
    `find_zero` to find the roots, when possible, **if** the interval is specified as an `Interval` object,
    as created by `-1..1`, say.


For example, this function (due to `@truculentmath`) is particularly tricky, as it is positive at every floating point number, but has two zeros (the vertical asymptote at `15//11` is only negative within adjacent floating point values):

```
julia> using IntervalArithmetic, IntervalRootFinding, Roots

julia> g(x) = x^2 + 1 +log(abs( 11*x-15 ))/99
g (generic function with 1 method)

julia> find_zeros(g, -3, 3)
Float64[]

julia> IntervalRootFinding.roots(g, -3..3, IntervalRootFinding.Bisection)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.36363, 1.36364], :unknown)
```

A less extreme usage might be the following, where `unique` indicates Bisection could be useful and indeed `find_zeros` will identify these values:

```
julia> g(x) = exp(x) - x^5
g (generic function with 1 method)

julia> rts = IntervalRootFinding.roots(g, -20..20)
2-element Vector{Root{Interval{Float64}}}:
 Root([12.7132, 12.7133], :unique)
 Root([1.29585, 1.29586], :unique)

julia> find_zeros(g, -20, 20)
2-element Vector{Float64}:
  1.2958555090953687
 12.713206788867632
```

"""
function find_zeros(f, a, b=nothing; no_pts=12, k=8, naive=false, kwargs...)
    if b === nothing
        a0, b0 = map(float, _extrema(a))
    else
        a0, b0 = promote(float(a), float(b))
    end
    a0 = isinf(a0) ? nextfloat(a0) : a0
    b0 = isinf(b0) ? prevfloat(b0) : b0

    # set tolerances if not specified
    fa0, fb0 = promote(float(f(a0)), float(f(b0)))
    d        = Dict(kwargs...)
    T, S     = real(eltype(a0)), real(eltype(fa0))
    xatol::T = get(d, :xatol, eps(one(T))^(4 / 5) * oneunit(T))
    xrtol    = get(d, :xrtol, eps(one(T)) * one(T))
    atol::S  = get(d, :atol, eps(float(S)) * oneunit(S))
    rtol     = get(d, :rtol, eps(float(S)) * one(S))

    zs = T[]  # collect zeros

    # check endpoints for exact zeros, then narrow
    abs(fa0) * oneunit(T) / oneunit(S) <= 8 * eps(a0) && push!(zs, a0)
    abs(fb0) * oneunit(T) / oneunit(S) <= 8 * eps(b0) && push!(zs, b0)
    a0 = find_non_zero(f, a0, b0, xatol, xrtol, atol, rtol)
    b0 = find_non_zero(f, b0, a0, xatol, xrtol, atol, rtol)

    (isnan(a0) || isnan(b0)) && throw(DomainError("no non-zero initial points found."))
    _fz!(zs, f, a0, b0, no_pts, k)  # initial zeros

    ints = Interval{T}[] # collect subintervals
    !naive &&
        !isempty(zs) &&
        make_intervals!(ints, f, a0, b0, zs, 1, xatol, xrtol, atol, rtol)

    nzs = T[]
    cnt = 0

    while !naive && !isempty(ints)
        cnt += 1
        i = pop!(ints)
        # this is fussy. Ideally, we would want to explore intervals
        # with fewer points (each interval is already probed for k=4
        # bisections) but how many fewer?  We already had ~ (i.b-i.a) *
        # no_pts / (b-a) points. Would want an increased density so n
        # > (i.b - i.a) / (b - a) * no_pts but this doesn't perform as
        # well as we might expect

        # sub_no_pts = ceil(Int, (i.b - i.a) / (b-a) * no_pts * 2^(i.depth))
        #sub_no_pts <= 2 && continue  # stop on depth, always divide if roots
        #sub_no_pts = max(3, floor(Int, no_pts  / (2.0)^(i.depth)))
        sub_no_pts = floor(Int, no_pts / (2.0)^(i.depth))

        empty!(nzs)
        if sub_no_pts >= 2
            _fz!(nzs, f, i.a, i.b, sub_no_pts, k)
        end

        if !isempty(nzs)
            azs = filter(rt -> not_near(rt, zs, xatol, xrtol), nzs) # trim out nearby roots
            length(azs) == 0 && continue
            append!(zs, azs)
            sort!(zs)
            i.depth > 4 && continue
            make_intervals!(ints, f, i.a, i.b, azs, i.depth + 1, xatol, xrtol, atol, rtol)
        end
    end

    length(zs) <= 1 && return zs
    sort!(zs)

    # may identify same zero with nearby values along the way
    # this trims out with a relaxed tolerance on how close
    # nearby roots can be. Default is epsilon^(2/5) ≈ 5e-7
    inds = Int[1]
    z1 = first(zs)
    for i in 2:length(zs)
        z2 = zs[i]
        if !approx_close(z1, z2, xatol, xrtol)
            push!(inds, i)
        end
        z1 = z2
    end

    return zs[inds]
end

# solve interface
"""
    AllZeros

Type to indicate to `solve` that `find_zeros` should be used to solve the given  `ZeroProblem`.

## Example

```
julia> Z = ZeroProblem(cos, (0, 2pi));

julia> solve(Z, AllZeros())
2-element Vector{Float64}:
 1.5707963267948966
 4.71238898038469
```

"""
struct AllZeros <: AbstractUnivariateZeroMethod end
function solve(𝑭𝑿::ZeroProblem, ::AllZeros; kwargs...)
    F, x₀ = 𝑭𝑿.F, 𝑭𝑿.x₀
    find_zeros(F, x₀; kwargs...)
end
